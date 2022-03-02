#!/usr/bin/env python


# MIT License
#
# Copyright (c) 2021-2022 Moritz E. Beber
# Copyright (c) 2021 Harshil Patel
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


"""Provide a tool for fetching experiment metadata from mixed database accessions."""


from __future__ import annotations

import argparse
import asyncio
import logging
import re
from csv import QUOTE_NONNUMERIC
from dataclasses import dataclass
from functools import partial
from io import StringIO
from itertools import chain
from pathlib import Path
from typing import (
    Iterable,
    Dict,
    Union,
    Optional,
    List,
    Pattern,
    ClassVar,
    Callable,
    Any,
)

import aiometer
import httpx
import pandas as pd
import rich.progress as rprog
import tenacity
from rich.logging import RichHandler


logger = logging.getLogger()


# Example ids supported by this script
SRA_IDS = (
    "PRJNA63463",
    "SAMN00765663",
    "SRA023522",
    "SRP003255",
    "SRR390278",
    "SRS282569",
    "SRX111814",
)
ENA_IDS = (
    "PRJEB7743",
    "SAMEA3121481",
    "ERA2421642",
    "ERP120836",
    "ERR674736",
    "ERS4399631",
    "ERX629702",
)
DDBJ_IDS = (
    "PRJDB4176",
    "SAMD00114846",
    "DRA008156",
    "DRP004793",
    "DRR171822",
    "DRS090921",
    "DRX162434",
)
GEO_IDS = ("GSE18729", "GSM465244")


@dataclass(frozen=True)
class SupportedDatabaseAccessions:
    """Define a value object for accessions sorted by type."""

    ena: List[str]
    sra: List[str]
    geo: List[str]

    _SRA_ID: ClassVar[Pattern] = re.compile(r"^(PRJNA|SAMN|SR[APRSX])\d+$")
    _ENA_ID: ClassVar[Pattern] = re.compile(r"^(PRJEB|SAMEA|ER[APRSX])\d+$")
    _DDBJ_ID: ClassVar[Pattern] = re.compile(r"^(PRJDB|SAMD|DR[APRSX])\d+$")
    _GEO_ID: ClassVar[Pattern] = re.compile(r"^G(PL|SM|SE|DS)\d+$")

    @classmethod
    def make_from_mixed(cls, accessions: Iterable[str]) -> SupportedDatabaseAccessions:
        """Create an instance from a mixed bunch of accessions."""
        ena = []
        sra = []
        geo = []
        for acc in accessions:
            if cls._ENA_ID.match(acc):
                ena.append(acc)
            elif cls._SRA_ID.match(acc) or cls._DDBJ_ID.match(acc):
                sra.append(acc)
            elif cls._GEO_ID.match(acc):
                geo.append(acc)
            else:
                raise ValueError(f"Unrecognized accession: {acc}")
        return SupportedDatabaseAccessions(ena=ena, sra=sra, geo=geo)


class RemoteFetcher:
    """Define a service for fetching remote resources concurrently and resiliently."""

    @classmethod
    async def fetch_all(
        cls,
        accessions: Iterable[str],
        url: Union[httpx.URL, str],
        query_kwargs: Dict[str, Union[int, str]],
        description: str,
        key: str,
        callback: Callable,
        **kwargs,
    ) -> List[Any]:
        """Fetch resources in a concurrent manner."""
        async with httpx.AsyncClient(
            base_url=url, params=query_kwargs
        ) as client:  # type: httpx.AsyncClient
            requests = cls._build_requests(client, accessions, key)
            result = await cls._fetch_concurrently(
                client, requests, description, callback=callback, **kwargs
            )
        return result

    @classmethod
    def _build_requests(
        cls,
        client: httpx.AsyncClient,
        accessions: Iterable[str],
        key: str,
    ) -> List[httpx.Request]:
        """Create one request object per accession."""
        return [
            client.build_request(method="GET", url="", params={key: acc})
            for acc in accessions
        ]

    @classmethod
    def _progress_bar(cls) -> rprog.Progress:
        """Create a rich progress bar."""
        return rprog.Progress(
            rprog.TextColumn("[bold blue]{task.description}", justify="right"),
            rprog.BarColumn(bar_width=None),
            "[progress.percentage]{task.completed}/{task.total}({task.percentage:>3.1f}%)",
            " ",
            rprog.TimeRemainingColumn(),
            " ",
            rprog.TimeElapsedColumn(),
        )

    @classmethod
    async def _fetch_concurrently(
        cls,
        client: httpx.AsyncClient,
        requests: List[httpx.Request],
        description: str,
        callback: Callable,
        max_concurrency: int = 10,
        max_per_second: int = 10,
        **kwargs,
    ) -> List[Any]:
        """Perform all given requests concurrently."""
        result = []
        with cls._progress_bar() as pbar:
            task = pbar.add_task(description=description, total=len(requests))
            async with aiometer.amap(
                partial(cls._fetch, client),
                requests,
                max_at_once=max_concurrency,
                max_per_second=max_per_second,
            ) as results:
                async for data in results:
                    result.append(callback(data, **kwargs))
                    pbar.update(task, advance=1)
        return result

    @staticmethod
    @tenacity.retry(
        wait=tenacity.wait_random_exponential(),
        stop=tenacity.stop_after_attempt(5),
        reraise=True,
        retry=tenacity.retry_if_exception_type(httpx.HTTPError),
        before=tenacity.before_log(logger, logging.DEBUG),
    )
    async def _fetch(
        client: httpx.AsyncClient,
        request: httpx.Request,
    ) -> StringIO:
        """Perform a request in an asynchronous and resilient manner."""
        response = await client.send(request, stream=True)
        response.raise_for_status()
        if response.status_code == 204:
            if "accession" in request.url.params:
                accession = request.url.params["accession"]
            elif "term" in request.url.params:
                accession = request.url.params["term"]
            else:
                accession = "?"
            raise ValueError(
                f"There is no content for accession {accession}. "
                f"Maybe you lack the right permissions?"
            )
        buffer = StringIO()
        async for chunk in response.aiter_text():
            buffer.write(chunk)
        buffer.seek(0)
        return buffer


class ExperimentAccessionService:
    """Define a service for fetching ENA experiment metadata from mixed accessions."""

    _GSM_PATTERN = re.compile(r"\^SAMPLE = (GSM\d+)")
    # List of metadata fields for ENA runs.
    # Full list of accepted fields can be obtained here:
    # https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&format=tsv&result=read_run
    ENA_METADATA_FIELDS = (
        "accession",
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "study_accession",
        "secondary_study_accession",
        "parent_study",
        "submission_accession",
        "run_alias",
        "experiment_alias",
        "sample_alias",
        "study_alias",
        "library_layout",
        "library_selection",
        "library_source",
        "library_strategy",
        "library_name",
        "instrument_model",
        "instrument_platform",
        "base_count",
        "read_count",
        "tax_id",
        "scientific_name",
        "sample_title",
        "experiment_title",
        "study_title",
        "description",
        "sample_description",
        "fastq_md5",
        "fastq_bytes",
        "fastq_ftp",
        "fastq_galaxy",
        "fastq_aspera",
    )

    @classmethod
    async def fetch_runinfo(cls, accessions: Iterable[str]) -> pd.DataFrame:
        """Return experiment metadata from the given accessions."""
        logger.info("Sort given accessions.")
        supported = SupportedDatabaseAccessions.make_from_mixed(accessions)
        logger.info("Assemble unique experiment accesssions.")
        experiments = await cls._unique_experiments(supported)
        logger.info("Fetch ENA experiment metadata.")
        return (
            (
                await cls._accession2erx(
                    experiments,
                    fields=cls.ENA_METADATA_FIELDS,
                    description="ENA Metadata",
                )
            )
            .drop_duplicates(subset=["experiment_accession", "run_accession"])
            .sort_values(
                by=["experiment_accession", "run_accession"], ignore_index=True
            )
        )

    @classmethod
    async def _unique_experiments(cls, accessions: SupportedDatabaseAccessions):
        """Return all unique experiment accessions from the given ones."""
        experiments = []
        if accessions.ena:
            df = await cls._accession2erx(accessions.ena)
            logger.debug(df.to_string())
            experiments.append(df["experiment_accession"])
        if accessions.sra:
            df = await cls._accession2srx(accessions.sra)
            logger.debug(df.to_string())
            experiments.append(df["Experiment"])
        if accessions.geo:
            df = await cls._geo2srx(accessions.geo)
            logger.debug(df.to_string())
            experiments.append(df["Experiment"])
        return set(chain.from_iterable(experiments))

    @classmethod
    async def _accession2srx(
        cls,
        accessions: Iterable[str],
        url: Union[
            httpx.URL, str
        ] = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi",
        query_kwargs: Optional[Dict[str, Union[int, str]]] = None,
        description: str = "➜ SRX",
    ) -> pd.DataFrame:
        """Resolve SRA, DDBJ, and GSM accessions to SRA experiments."""
        if query_kwargs is None:
            query_kwargs = {"save": "efetch", "db": "sra", "rettype": "runinfo"}
        tables = await RemoteFetcher.fetch_all(
            accessions=accessions,
            url=url,
            query_kwargs=query_kwargs,
            description=description,
            key="term",
            callback=cls._parse_table,
            delimiter=",",
        )
        return pd.concat(tables, ignore_index=True)

    @classmethod
    async def _geo2srx(
        cls,
        accessions: Iterable[str],
        url: Union[httpx.URL, str] = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
        query_kwargs: Optional[Dict[str, Union[int, str]]] = None,
        description: str = "➜ GSM",
    ) -> pd.DataFrame:
        """Resolve GEO accessions to SRA experiments."""
        if query_kwargs is None:
            query_kwargs = {"targ": "gsm", "view": "data", "form": "text"}
        gsm = await RemoteFetcher.fetch_all(
            accessions=accessions,
            url=url,
            query_kwargs=query_kwargs,
            description=description,
            key="acc",
            callback=cls._parse_gsm,
        )
        return await cls._accession2srx(gsm)

    @classmethod
    async def _accession2erx(
        cls,
        accessions: Iterable[str],
        fields: Iterable[str] = ("experiment_accession", "run_accession"),
        url: Union[httpx.URL, str] = "https://www.ebi.ac.uk/ena/portal/api/filereport",
        query_kwargs: Optional[Dict[str, Union[int, str]]] = None,
        description: str = "➜ ERX",
    ) -> pd.DataFrame:
        """Resolve ENA accessions to ENA experiments."""
        if query_kwargs is None:
            query_kwargs = {"result": "read_run", "limit": 0}
        query_kwargs["fields"] = ",".join(fields)
        tables = await RemoteFetcher.fetch_all(
            accessions=accessions,
            url=url,
            query_kwargs=query_kwargs,
            description=description,
            key="accession",
            callback=cls._parse_table,
            delimiter="\t",
        )
        return pd.concat(tables, ignore_index=True)

    @classmethod
    def _parse_table(cls, data: StringIO, delimiter: str) -> pd.DataFrame:
        """Parse a table from the given text."""
        return pd.read_csv(data, sep=delimiter, header=0, index_col=False)

    @classmethod
    def _parse_gsm(cls, data: StringIO) -> List[str]:
        """Parse all GSM accessions from the given text."""
        return cls._GSM_PATTERN.findall(data.getvalue())


@tenacity.retry(
    wait=tenacity.wait_random_exponential(),
    stop=tenacity.stop_after_attempt(5),
    reraise=True,
    retry=tenacity.retry_if_exception_type(httpx.HTTPError),
    before=tenacity.before_log(logger, logging.DEBUG),
)
def fetch_metadata_fields(
    url: Union[httpx.URL, str] = "https://www.ebi.ac.uk/ena/portal/api/returnFields",
    query_kwargs: Optional[Dict[str, Union[int, str]]] = None,
) -> Iterable[str]:
    """Fetch the allowed metadata fields from the ENA API."""
    if query_kwargs is None:
        query_kwargs = {"dataPortal": "ena", "format": "tsv", "result": "read_run"}
    response = httpx.get(url=url, params=query_kwargs)
    response.raise_for_status()
    result = pd.read_csv(StringIO(response.text), sep="\t", header=0, index_col=False)
    return result["columnId"]


def validate_metadata_fields(fields: List[str]) -> None:
    """Expect that the given metadata fields are a subset of the allowed fields."""
    allowed = fetch_metadata_fields()
    if not set(fields).issubset(allowed):
        raise ValueError(
            f"Please provide valid ENA metadata fields.\n"
            f"Provided fields: {', '.join(fields)}\n"
            f"Allowed fields: {', '.join(allowed)}\n"
        )


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Download and create a run information metadata file from SRA / "
        "ENA / DDBJ / GEO identifiers.",
        # epilog="Example usage: python fetch_sra_runinfo.py <FILE_IN> <FILE_OUT>",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="File containing database identifiers, one per line.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Output file in tab-delimited format.",
    )
    parser.add_argument(
        "-m",
        "--ena-metadata-fields",
        type=str,
        default=",".join(ExperimentAccessionService.ENA_METADATA_FIELDS),
        help="Comma-separated list of ENA metadata fields to fetch.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(args)


def main(args: Optional[List[str]] = None) -> None:
    """Manage input and output as well as program flow."""
    args = parse_args(args)
    logging.basicConfig(
        level=args.log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(markup=True, rich_tracebacks=True)],
    )
    if not args.file_in.is_file():
        raise FileNotFoundError(f"The given input file {args.file_in} was not found!")
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    ExperimentAccessionService.ENA_METADATA_FIELDS = [
        field
        for token in args.ena_metadata_fields.split(",")
        if (field := token.strip())
    ]
    validate_metadata_fields(ExperimentAccessionService.ENA_METADATA_FIELDS)
    with args.file_in.open() as handle:
        accessions = [acc for line in handle.readlines() if (acc := line.strip())]
    metadata: pd.DataFrame = asyncio.run(
        ExperimentAccessionService.fetch_runinfo(accessions)
    )
    metadata.to_csv(
        args.file_out, sep="\t", header=True, index=False, quoting=QUOTE_NONNUMERIC
    )


if __name__ == "__main__":
    main()
