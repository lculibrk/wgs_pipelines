from genomicsapi import connection, q
import genomicsapi.dbs
import genomicsapi.translate

def describe(db, table):
    cursor = connection.cursor(buffered = True)
    genomicsapi.dbs.chdb(db, cursor)
    query = f"DESCRIBE {table}"
    cursor.execute(query)
    result = cursor.fetchall()
    return(result)

def simple_select(db, table, filter_column, filter_value):
    ## Initialize
    cursor = connection.cursor(buffered = True)
    genomicsapi.dbs.chdb(db, cursor)
    ## Sanitize filter_value
    filter_value = q(filter_value)
    ## Build the query
    query = f"SELECT * FROM {table} WHERE {filter_column} = {filter_value}"
    cursor.execute(query)
    result = cursor.fetchall()
    return(result)

def multi_select(db, table, filters):
    cursor = connection.cursor(buffered = True)
    genomicsapi.dbs.chdb(db, cursor)
    if len(filters.keys()) == 1:
        return(simple_select(db, table, list(filters.keys())[0], list(filters.values())[0]))
    columns = list(filters.keys())
    filters = [f"{k} = {q(v)}" for k,v in filters.items()]
    filterstring = " AND ".join(filters)
    query = f"SELECT * FROM {table} WHERE {filterstring}"
    print(f"Trying {query}")
    cursor.execute(query)
    result = cursor.fetchall()
    return(result)

def parent_ids(in_id, db = "genomicsdb"):
    ## Prefix determines the database table for the ID
    prefix = in_id[:3]
    if prefix == "S":
        idtype = "samples"
    elif prefix == "R":
        idtype = "runs"
    elif prefix == "P":
        idtype = "studies"
    #print(id_dict)
    #print(prefix)
    #print(idtype)
    id_dict = {"studies":None, "samples":None, "runs":None}
    id_dict[idtype] = in_id
    #print(id_dict)
    ## End recursiveness
    if idtype == "studies":
        return(id_dict)
    ## Numeric ID
    num = int(in_id[3:])
    result = genomicsapi.select.simple_select(db, idtype, "id", num)
    if not result:
        raise ValueError("Supplied ID {in_id} has no match in database")
    result = result[0]
    table_descr = genomicsapi.select.describe(db, idtype)
    table_descr = [element[0] for element in table_descr]
    result = dict(zip(table_descr, result))
    id_dict["studies"] = genomicsapi.translate.idtostring(result["study_id"], "P")
    if idtype == "runs":
        id_dict["samples"] = genomicsapi.translate.idtostring(result["sample_id"], "S")
    return(id_dict)