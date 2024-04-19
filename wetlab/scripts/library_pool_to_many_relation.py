import wetlab.models


def run(f_name):
    """This script is part of the issue "#180,when deleting run , pool and
        library_preparations are also deleted" on Class LibraryPool.
        The first part of the script fetch the existing data on run_process_id
        and create a file containing the pk of the libraryPool instance and
        the pk of the run_process_id
        The second part of the script fetch the data from the file and add the
        run_process_pk to the run_process field of the LibraryPool instance
    """
    if hasattr(wetlab.models.LibraryPool, "run_process_id"):
        with open (f_name, "w") as fo:
            for library_pool in wetlab.models.LibraryPool.objects.all():
                # if library_pool.run_process_id:
                #    library_pool.run_process.add(library_pool.run_process_id)
                if library_pool.run_process_id is not None:
                    fo.write(str(library_pool.id) + "," + str(library_pool.run_process_id.id) + "\n")
    if hasattr(wetlab.models.LibraryPool, "run_process"):
        with open (f_name, "r") as fh:
            for line in fh:
                library_pool_pk, run_process_pk = line.strip().split(",")
                try:
                    library_pool = wetlab.models.LibraryPool.objects.get(id=library_pool_pk)
                    run_process = wetlab.models.RunProcess.objects.get(id=run_process_pk)
                except Exception as e:
                    print(f"Error: {e}")
                    continue
                if not library_pool.run_process.filter(id=run_process_pk).exists():
                    library_pool.run_process.add(run_process)

            
    
