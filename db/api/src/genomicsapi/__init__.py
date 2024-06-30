import os
import mysql.connector

def q(input_string):
    if isinstance(input_string, str):
        if not ((input_string.startswith('"') and input_string.endswith('"')) or
                (input_string.startswith("'") and input_string.endswith("'"))):
            return f'"{input_string}"' if "'" in input_string else f"'{input_string}'"
    elif input_string is None:
        return("NULL")
    return input_string

connection = mysql.connector.connect(host = os.environ["DBHOST"], user = os.environ["DBUSER"], passwd = os.environ["DBPASS"])