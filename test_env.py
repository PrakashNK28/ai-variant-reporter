import os
from dotenv import load_dotenv

# Force load .env from current folder
load_dotenv(dotenv_path=".env")

key = os.getenv("ANTHROPIC_API_KEY")

print("API KEY:", key)