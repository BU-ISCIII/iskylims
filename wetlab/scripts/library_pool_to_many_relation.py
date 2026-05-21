import wetlab.models


def run(f_name):
    """Upgrade: 3.0.0 -> 3.1.0
    This script is part of the issue "#180,when deleting run , pool and
    library_preparations are also deleted" on Class LibraryPool.
    This script expects a file created externally that contains the pk of the
    LibraryPool instance and the pk of the run_process_id.
    Export example:
        mysql -u <db_user> -p -h <db_host> -D <db_name> \
          -e "SELECT id, run_process_id_id FROM wetlab_library_pool" \
          > /tmp/library_pool_run_process.tsv
    The script fetches the data from the file and adds the
    run_process_pk to the run_process field of the LibraryPool instance
    """
    if hasattr(wetlab.models.LibraryPool, "run_process"):
        with open(f_name, "r") as fh:
            lines = fh.readlines()
        heading = lines[0].strip().split("\t")
        if "run_process_id_id" in heading:
            r_proc_idx = heading.index("run_process_id_id")
            pk_idx = 0
            for line in lines[1:]:
                l_data = line.strip().split("\t")
                if l_data[r_proc_idx] == "NULL":
                    continue
                library_pool_pk = l_data[pk_idx]
                run_process_pk = l_data[r_proc_idx]

                try:
                    library_pool = wetlab.models.LibraryPool.objects.get(
                        id=library_pool_pk
                    )
                    run_process = wetlab.models.RunProcess.objects.get(
                        id=run_process_pk
                    )
                except Exception as e:
                    print(f"Error: {e}")
                    continue
                if not library_pool.run_process.filter(id=run_process_pk).exists():
                    library_pool.run_process.add(run_process)
                    print(
                        f"LibraryPool: {library_pool_pk} added run_process: {run_process_pk}"
                    )
        else:
            print("Data migration for LibraryPool was already done")
    else:
        print("LibraryPool model does not have run_process field")
