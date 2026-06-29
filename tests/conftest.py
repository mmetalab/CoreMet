"""Shared fixtures for CoreMet test suite."""

import pytest
import sys
from pathlib import Path

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture(scope="session")
def app():
    """Create a Dash app instance for testing."""
    from app.main import create_app
    application = create_app()
    return application


@pytest.fixture(scope="session")
def client(app):
    """Flask test client."""
    return app.server.test_client()
