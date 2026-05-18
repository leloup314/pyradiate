import datetime
import pytest
from pathlib import Path
from pyradiate.tools import ensdf_lib_updater

ENSDF_ARCHIVE_TEST_DATE = datetime.date(2026, 5, 1)


def test_ensdf_lib_updater_retrieves_latest_date():
    latest_archive = ensdf_lib_updater.get_latest_archive_date()
    assert isinstance(latest_archive, datetime.date)
    assert latest_archive >= ENSDF_ARCHIVE_TEST_DATE


def test_endsdf_library_updater(tmpdir):
    ensdf_file = tmpdir / Path("ensdf.zip")
    ensdf_lib_updater.download_ensdf_archive(
        archive_date=ENSDF_ARCHIVE_TEST_DATE,
        archive_file=ensdf_file,
    )
    assert ensdf_file.isfile()
