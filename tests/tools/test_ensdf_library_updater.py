import datetime
import pytest
import random
import zipfile
from pathlib import Path
from pyradiate.tools import ensdf_lib_updater

ENSDF_ARCHIVE_TEST_DATE = datetime.date(2026, 5, 1)


@pytest.fixture(scope="session")
def session_tmp_path(tmp_path_factory):
    return tmp_path_factory.mktemp("data")


@pytest.fixture()
def ensdf_file(session_tmp_path):
    return session_tmp_path / "ensdf.zip"


@pytest.fixture()
def archive_path(session_tmp_path):
    return session_tmp_path / "archive_path"


def test_get_latest_archive_date():
    latest_archive = ensdf_lib_updater.get_latest_archive_date()
    assert isinstance(latest_archive, datetime.date)
    assert latest_archive >= ENSDF_ARCHIVE_TEST_DATE


@pytest.mark.dependency()
def test_downloading_archive(ensdf_file):
    assert not ensdf_file.is_file()
    ensdf_lib_updater.download_ensdf_archive(
        archive_date=ENSDF_ARCHIVE_TEST_DATE,
        archive_file=ensdf_file,
    )
    assert ensdf_file.is_file()


@pytest.mark.dependency(depends=["test_downloading_archive"])
def test_unpacking_archive(ensdf_file, archive_path):
    assert not archive_path.exists()
    ensdf_lib_updater.unpack_endsf_archive(archive_file=ensdf_file, archive_path=archive_path)
    assert archive_path.exists()
    expected_files = set(f"ensdf.{i:03d}" for i in range(1, 301))
    for _file in archive_path.iterdir():
        assert _file.name in expected_files


@pytest.mark.dependency(depends=["test_unpacking_archive"])
def test_unpacking_archive_fails_if_file_missing(ensdf_file, archive_path, session_tmp_path):
    bad_ensdf_file = session_tmp_path / Path("ensdf_bad.zip")
    with zipfile.ZipFile(bad_ensdf_file, "w") as bad:
        with zipfile.ZipFile(ensdf_file, "r") as good:
            content = good.infolist()
            bad_idx = random.randint(0, len(content))
            for idx, item in enumerate(content):
                if idx == bad_idx:
                    continue
                bad.writestr(item, good.read(item.filename))

    archive_path = session_tmp_path / Path("bad_archive_path")
    with pytest.raises(RuntimeWarning):
        ensdf_lib_updater.unpack_endsf_archive(archive_file=bad_ensdf_file, archive_path=archive_path)
    assert not archive_path.exists()
