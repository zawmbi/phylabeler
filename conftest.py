"""Pytest glue for the hand-rolled test scripts.

test_monophyly.py and test_archive_io.py are written to run standalone -
`python test_monophyly.py` - with a main() that builds one shared object and
passes it to each test in turn. Pytest collects those same `test_*` functions
but reads their parameters as fixture requests, so every test that wanted the
prepared taxonomy or fixture directory errored out with "fixture 'db' not
found" before running a line.

The fixtures below hand pytest the very objects main() builds, so both entry
points exercise the same code with the same inputs and neither one has to
change.
"""

import pytest


@pytest.fixture(scope="session")
def db():
    """The synthetic taxonomy from test_monophyly.main()."""
    import test_monophyly
    return test_monophyly.build_db()


@pytest.fixture
def path(tmp_path):
    """Fixture directory from test_archive_io.main(), one per test.

    Function-scoped on purpose: these tests write caches and extract archives
    into this directory, so sharing one across the module would let an earlier
    test decide what a later one sees. tmp_path stands in for the script's
    TemporaryDirectory and pytest cleans it up.
    """
    import test_archive_io
    return test_archive_io.build_fixtures(str(tmp_path))


@pytest.hookimpl(hookwrapper=True)
def pytest_pyfunc_call(pyfuncitem):
    """Turn the scripts' check() tally into a real pytest failure.

    check() records a mismatch and prints it, but never raises - that is what
    lets the standalone runner report every failure in one go instead of
    stopping at the first. Under pytest the same design would report a green
    test with output full of FAILs, which is worse than the errors it
    replaced. Comparing the module's counter either side of the call keeps
    run-everything behaviour and still fails the test if anything missed.

    This is a call-phase wrapper rather than an autouse fixture on purpose:
    pytest.fail() during fixture teardown is reported as an ERROR, which is
    exactly the category this file exists to empty. Raised here it lands as
    an honest FAILURE.
    """
    module = pyfuncitem.module
    before = getattr(module, "_failed", None)
    outcome = yield
    outcome.get_result()      # let a genuine exception win over the tally
    after = getattr(module, "_failed", None)
    if before is not None and after is not None and after > before:
        raise AssertionError(
            f"{after - before} check(s) failed in {pyfuncitem.name} - see the "
            f"captured output above")
