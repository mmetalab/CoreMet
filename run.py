"""
Application entry point for CoreMet Web Server
"""

import os
import sys
from pathlib import Path

# Add the project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Load .env file (must be before any service imports)
from dotenv import load_dotenv
load_dotenv(project_root / '.env')

from app.main import create_app

# Create application instance
app = create_app()

# Server object for WSGI deployment
server = app.server

if __name__ == '__main__':
    # Get configuration from environment
    config_name = os.getenv('FLASK_ENV', 'development')
    
    # Run the application
    app.run(
        debug=(config_name == 'development'),
        host='0.0.0.0',
        port=int(os.getenv('PORT', 8080)),
        use_reloader=False,  # Reloader causes hang with eager page imports
    )
