def chdb(dbname, cursor):
    cursor.execute(f"USE {dbname};")
    return(None)