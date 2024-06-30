from genomicsapi import connection, q
import genomicsapi.dbs
import genomicsapi.select

def generic_insert(insert_keys, table, db = "genomicsdb_devel"):
    ## Initialize
    cursor = connection.cursor(buffered = True)
    genomicsapi.dbs.chdb(db, cursor)
    ## Check types
    if not isinstance(insert_keys, dict):
        raise(TypeError("insert_keys must be a dict"))
    if "rname" not in insert_keys.keys():
        raise(ValueError("rname must be a key in insert_keys"))
    ## unpack rname
    rname = insert_keys["rname"]
    ## Check if what we're inserting already exists
    result = genomicsapi.select.multi_select(db = db, table = table, filters = insert_keys)
    ## If not, then do the insert
    if result:
        print(f"Value {rname} in column rname already exists. Skipping this add.")
    else:
        ## Create the colnames/values strings
        ## Need to keep things sorted so we do it this way
        ## Convert everything to a string, things that were strings before get quoted to play nice with mySQL
        colkeys = list(insert_keys.keys())
        colvals = [str(q(insert_keys[k])) for k in colkeys]
        colkeys = ", ".join(colkeys)
        colvals = ", ".join(colvals)
        ## Create the query and add the line
        query = f"INSERT INTO {table}({colkeys}) VALUES ({colvals})"
        print(f"Performing operation:\n{query};")
        cursor.execute(query)
        connection.commit()
        ## Now get the result of what we just inserted
        result = genomicsapi.select.multi_select(db = db, table = table, filters = insert_keys)
    return(result)

def analysis_insert(insert_keys, table, db = "genomicsdb_devel"):
    ## Initialize
    cursor = connection.cursor(buffered = True)
    genomicsapi.dbs.chdb(db, cursor)
    ## Check types
    if not isinstance(insert_keys, dict):
        raise(TypeError("insert_keys must be a dict"))
    ## Determine the ID we're dealing with
    if "runs_id" in insert_keys:
        tbl = "runs_id"
        tbl_id = insert_keys["runs_id"]
    elif "samples_id" in insert_keys:
        tbl = "samples_id"
        tbl_id = insert_keys["samples_id"]
    elif "studies_id" in insert_keys:
        tbl = "studies_id"
        tbl_id = insert_keys["studies_id"]
    ## Check if what we're inserting already exists
    result = genomicsapi.select.multi_select(db = db, table = table, filters = {"pipeline_name":insert_keys["pipeline_name"], "pipeline_version":insert_keys["pipeline_version"], tbl:tbl_id})
    ## If not, then do the insert
    if result:
        print(f"Row with pipeline name {insert_keys['pipeline_name']}, version {insert_keys['pipeline_version']} for {tbl} {tbl_id} already exists. Skipping this add.")
    else:
        ## Create the colnames/values strings
        ## Need to keep things sorted so we do it this way
        ## Convert everything to a string, things that were strings before get quoted to play nice with mySQL
        colkeys = list(insert_keys.keys())
        colvals = [str(q(insert_keys[k])) for k in colkeys]
        colkeys = ", ".join(colkeys)
        colvals = ", ".join(colvals)
        ## Create the query and add the line
        query = f"INSERT INTO {table}({colkeys}) VALUES ({colvals})"
        print(f"Performing operation:\n{query};")
        cursor.execute(query)
        connection.commit()
        ## Now get the result of what we just inserted
        result = genomicsapi.select.multi_select(db = db, table = table, filters = {"pipeline_name":insert_keys["pipeline_name"], "pipeline_version":insert_keys["pipeline_version"], tbl:tbl_id})
    return(result)
