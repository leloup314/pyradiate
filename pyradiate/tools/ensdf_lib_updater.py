import datetime
import zipfile
import yaml

from  urllib.request import urlopen
from html.parser import HTMLParser
from tqdm import tqdm

from pyradiate import logger, ensdf_file, ensdf_path, ensdf_config


# URL of archives for nuclide libraries and known existing date at the time of writing
NBL_NNDC_ENDSF_ARCHIVE_URL = "https://www.nndc.bnl.gov/ensdfarchivals/"
FALLBACK_LIB_DATE = datetime.date(year=2025, month=8, day=4)


class NBLDataParser(HTMLParser):
    data = None
    def handle_data(self, data):
        self.data = data


def get_latest_archive_date():
    """Returns date of latest ENSDF archive as datetime object"""
    
    latest_date = None
    parser = NBLDataParser()
    
    # Read entire page
    with urlopen(NBL_NNDC_ENDSF_ARCHIVE_URL) as page:
        content = [l.decode("utf-8").strip() for l in page.readlines()]
        for i, cntnt in enumerate(content):
            parser.feed(cntnt)
            # Next entry is date of latest zip archive
            # Date format is "YYYY/MM/DD"; needs to be brought in isofromat "YYYY-MM-DD"
            if parser.data == "Last modified:":
                parser.feed(content[i+1])
                latest_date = parser.data.replace('/', '-')
                break

    parser.close()
    
    try:
        return datetime.date.fromisoformat(latest_date)
    except ValueError as e:
        logger.warning(f"Could not determine date of latest ENSDF archive.")
        return None
    

def build_nbl_ensdf_archive_url(archive_date):
    """Build the URL under which the ENDSF archives are found for a specific *date*"""
    archive_file = f"ensdf_{archive_date.strftime('%y%m%d')}.zip"  # e.g. "ensdf_250804.zip"
    dist_folder = f"dist{archive_date.strftime('%y')}"  # e.g. "dist25" for year 2025
    return f"{NBL_NNDC_ENDSF_ARCHIVE_URL}/distributions/{dist_folder}/{archive_file}"


def download_endsf_archive(archive_date, archive_file=ensdf_file, chunk_size=512*1024):
    """Downloads an ENSDF archive from *date* to the pyradiate.data_path folder"""
    
    archive_url = build_nbl_ensdf_archive_url(archive_date)
    
    with urlopen(archive_url) as archive_stream:

        # Verify returncode
        assert archive_stream.status == 200, f"Incorrect return code {archive_stream.status}"

        # Extract archive info
        archive_info = archive_stream.headers
        archive_n_bytes = int(archive_info['Content-Length'])
        
        # Logging
        logger.debug(f"Retrieved archive information:\n{archive_info}")
        logger.info(f"Retrieving {archive_n_bytes/1024**2:.2f} MB ENSDF archive from {archive_date}...")

        with open(archive_file, 'wb') as archive_local:
            pbar = tqdm(total=archive_n_bytes, unit='bytes', unit_scale=True)
            while chunk := archive_stream.read(chunk_size): 
                archive_local.write(chunk)
                pbar.update(len(chunk))

        assert archive_file.stat().st_size == archive_n_bytes, "Size of local and remote ENSDF archives differ"

        # Write ENSDF archive config file
        archive_config = {
            "date": archive_info["Last-Modified"],
            "size": archive_n_bytes
        } 
        with open(ensdf_config, 'w') as conf:
            yaml.safe_dump(archive_config, conf)


        
def unpack_endsf_archive(archive_file=ensdf_file, archive_path=ensdf_path):
    """"""
    logger.debug(f"Unpacking ENSDF archive {archive_file} at {archive_path}")
    with zipfile.ZipFile(archive_file, 'r') as archive_zip:
        ensdf_contents = archive_zip.namelist()
        ensdf_missing = [i for i in range(1, 301) if f'ensdf.{i:03d}' not in ensdf_contents]
        if ensdf_missing:
            raise RuntimeWarning(f"ENSDF archive missing entries: {' ,'.join(f'ensdf.{i:03d}' for i in ensdf_missing)}")
        archive_zip.extractall(archive_path)
    logger.debug("ENSDF archive unpacked successfully")


def read_ensdf_archive_config(archive_config=ensdf_config):

    assert archive_config.exists(), "Archive file does not exists"
    
    with open(archive_config, 'r') as conf:
        res = yaml.safe_load(conf)
    
    res['date'] = datetime.datetime.strptime(res['date'], "%a, %d %b %Y %H:%M:%S %Z")
    
    return res

def main():

    # Latest ENSDF archive
    latest_ensdf = get_latest_archive_date()

    if latest_ensdf is None:
        logger.warning(f"Latest ENSDF archive date could not be read.")

    if not ensdf_file.is_file():
        logger.info("Required ENSDF archive file not found. Fetching archive...")

        if latest_ensdf is None:
            latest_ensdf = FALLBACK_LIB_DATE
            logger.warning(f"Using fallback ENSDF archive from {FALLBACK_LIB_DATE}")
        
        download_endsf_archive(archive_date=latest_ensdf, archive_file=ensdf_file)
        unpack_endsf_archive(archive_file=ensdf_file, archive_path=ensdf_path)
    
    else:
        archive_config = read_ensdf_archive_config(archive_config=ensdf_config)
        logger.info(f"Found ENSDF archive from {archive_config['date'].date()}")

        if latest_ensdf is not None and archive_config['date'].date() < latest_ensdf:
            logger.info(f"Updating ENSDF archive from {archive_config['date'].date()} with {latest_ensdf} version")
            download_endsf_archive(archive_date=latest_ensdf, archive_file=ensdf_file)
            unpack_endsf_archive(archive_file=ensdf_file, archive_path=ensdf_path)


if __name__ == "__main__":
    main()