import datetime
import zipfile
import yaml

from urllib.request import urlopen
from html.parser import HTMLParser
from tqdm import tqdm

from pyradiate import logger, ensdf_file, ensdf_path, ensdf_config
from pyradiate.tools.common import load_yaml, save_yaml


# URL of archives for nuclide libraries and known existing date at the time of writing
NBL_NNDC_ENDSF_ARCHIVE_URL = "https://www.nndc.bnl.gov/ensdfarchivals/"
FALLBACK_LIB_DATE = datetime.date(year=2026, month=5, day=1)


class NBLDataParser(HTMLParser):
    data = None

    def handle_data(self, data):
        self.data = data


def get_latest_archive_date() -> datetime.date | None:
    """Returns date of latest ENSDF archive as datetime object"""

    latest_date = None
    parser = NBLDataParser()

    # Read entire page
    with urlopen(NBL_NNDC_ENDSF_ARCHIVE_URL) as page:
        content = [line.decode("utf-8").strip() for line in page.readlines()]
        for i, cntnt in enumerate(content):
            parser.feed(cntnt)
            # Next entry is date of latest zip archive
            # Date format is "YYYY/MM/DD"; needs to be brought in isofromat "YYYY-MM-DD"
            # NOTE: the "Last modified" refers to the upload date. The"Last-Modified" field
            #       in archive_info refers to the actual file timestamp of the ZIP archive.
            if parser.data == "Last modified:":
                parser.feed(content[i + 1])
                latest_date = parser.data.replace("/", "-")
                break

    parser.close()

    try:
        return datetime.date.fromisoformat(latest_date)
    except ValueError as e:
        logger.warning(f"Could not determine date of latest ENSDF archive: {e!r}")
        return None


def build_nbl_ensdf_archive_url(archive_date) -> str:
    """Build the URL under which the ENDSF archives are found for a specific *date*"""
    archive_file = f"ensdf_{archive_date.strftime('%y%m%d')}.zip"  # e.g. "ensdf_250804.zip"
    dist_folder = f"dist{archive_date.strftime('%y')}"  # e.g. "dist25" for year 2025
    return f"{NBL_NNDC_ENDSF_ARCHIVE_URL}/distributions/{dist_folder}/{archive_file}"


def download_ensdf_archive(archive_date, archive_file=ensdf_file, chunk_size=512 * 1024) -> None:
    """Downloads an ENSDF archive from *date* to the pyradiate.data_path folder"""

    archive_url = build_nbl_ensdf_archive_url(archive_date)

    with urlopen(archive_url) as archive_stream:
        # Verify returncode
        assert archive_stream.status == 200, f"Incorrect return code {archive_stream.status}"

        # Extract archive info
        archive_info = archive_stream.headers
        archive_n_bytes = int(archive_info["Content-Length"])
        # Logging
        logger.debug(f"Retrieved archive information:\n{archive_info}")
        logger.info(f"Retrieving {archive_n_bytes / 1024**2:.2f} MB ENSDF archive from {archive_date}...")

        with open(archive_file, "wb") as archive_local:
            pbar = tqdm(total=archive_n_bytes, unit="bytes", unit_scale=True)
            while chunk := archive_stream.read(chunk_size):
                archive_local.write(chunk)
                pbar.update(len(chunk))

        assert archive_file.stat().st_size == archive_n_bytes, "Size of local and remote ENSDF archives differ"

    return archive_info


def unpack_endsf_archive(archive_file=ensdf_file, archive_path=ensdf_path) -> None:
    """"""
    logger.debug(f"Unpacking ENSDF archive {archive_file} at {archive_path}")
    with zipfile.ZipFile(archive_file, "r") as archive_zip:
        ensdf_contents = archive_zip.namelist()
        ensdf_missing = [i for i in range(1, 301) if f"ensdf.{i:03d}" not in ensdf_contents]
        if ensdf_missing:
            raise RuntimeWarning(f"ENSDF archive missing entries: {' ,'.join(f'ensdf.{i:03d}' for i in ensdf_missing)}")
        archive_zip.extractall(archive_path)
    logger.debug("ENSDF archive unpacked successfully")


def update_ensdf() -> None:
    """
    Function that
      - checks whether an ENSDF archive file exists and if not fetches the latest
      - updates the ENSDF archive file to the latest file if the local one is older
      - extracts / overwrites the local ENSDF files
      - upates the configuration accordingly
    All ENSDF-related files are located under *ensdf_path*
    """
    # Latest ENSDF archive
    latest_ensdf = get_latest_archive_date()

    if latest_ensdf is None:
        logger.warning(
            f"Latest ENSDF archive date could not be fetched. \
            Using fallback ENDSDF archive from {FALLBACK_LIB_DATE!r}"
        )
        latest_ensdf = FALLBACK_LIB_DATE

    update_ensdf_archive_file = False
    if not ensdf_file.is_file():
        update_ensdf_archive_file = True
        logger.info("Required ENSDF archive file not found. Fetching archive...")
        download_ensdf_archive(archive_date=latest_ensdf, archive_file=ensdf_file)
        unpack_endsf_archive(archive_file=ensdf_file, archive_path=ensdf_path)
    else:
        archive_config = load_yaml(ensdf_config)
        logger.info(f"Found ENSDF archive from {archive_config["date"]}")
        if datetime.date.fromisoformat(archive_config["date"]) < latest_ensdf:
            update_ensdf_archive_file = True
            logger.info(f"Updating ENSDF archive from {archive_config["date"]} with {latest_ensdf} version")
            download_ensdf_archive(archive_date=latest_ensdf, archive_file=ensdf_file)
            unpack_endsf_archive(archive_file=ensdf_file, archive_path=ensdf_path)

    # Write ENSDF archive config file
    if update_ensdf_archive_file:
        save_yaml(ensdf_config, {"date": latest_ensdf.isoformat()})


if __name__ == "__main__":
    update_ensdf()
