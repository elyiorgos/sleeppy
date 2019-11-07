from pytest import fixture
from importlib import resources


@fixture(scope="package")
def geneactiv_data():
    def _get_path(data):
        if data == "bin":
            with resources.path(
                "sleeppy.test.test_data", "test_geneactiv.bin"
            ) as file_path:
                return file_path
        elif data == "csv":
            with resources.path(
                "sleeppy.test.test_data", "test_geneactiv.csv"
            ) as file_path:
                return file_path

    return _get_path


@fixture(scope="package")
def activity_index_data():
    with resources.path(
        "sleeppy.test.test_data", "test_activity_index.h5"
    ) as file_path:
        return file_path


@fixture(scope="package")
def output_table():
    def _get_path(data):
        if data == "input":
            with resources.path(
                "sleeppy.test.test_data", "test_report_output_table.csv"
            ) as file_path:
                return file_path
        elif data == "expected":
            with resources.path(
                "sleeppy.test.test_data", "test_report_output_table_expected.csv"
            ) as file_path:
                return file_path

    return _get_path
