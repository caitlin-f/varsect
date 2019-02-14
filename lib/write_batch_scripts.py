"""
Writes the Slurm batch scripts for each stage of Varsect.
Script is called by 'Varsect_batch.py'
"""
def write_batch_files(args, file_sets):
    """ Write the slurm batch files """

    if args.M:
        __mapping_batch(args, file_sets)

    if args.D:
        __delly_batch(args, file_sets)

    if args.F:
        __freebayes_batch(args, file_sets)

    if args.G:
        __gatk_batch(args, file_sets)

    if args.P:
        __pindel_batch(args, file_sets)

    if args.C:
        __filter_batch(args, file_sets)
        # __intersect_batch(args, file_sets)
        # __collate_batch(args, file_sets)
        __matrix_batch(args, file_sets)

    if args.R:
        __raxml_batch(args, file_sets)

    if args.B:
        __mrbayes_batch(args, file_sets)

    if args.I:
        __samplot_batch(args, file_sets)

def __mapping_batch(args, file_sets):
    """ Writes batch script for mapping with bwa """
    with open("{}/1_Mapping/mapping_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name bwa\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task {}\n".format(args.t))
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%{}\n".format(len(file_sets.keys()), 48//args.t))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.map.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.map.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write('FWD=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/forward.reads )\n'.format(args.o))
        outfile.write('REV=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/reverse.reads )\n'.format(args.o))
        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("REF={}\n".format(args.r))
        outfile.write("THREADS={}\n\n".format(args.t))

        outfile.write('echo "Running bwa on ${SAMPLE}"\n\n')

        outfile.write("bash /home.roaming/s4097594/SV_pipeline/scripts/mapping.sh -t ${THREADS} -s ${SAMPLE} -r ${REF} -1 ${FWD} -2 ${REV} -o ${OUTDIR}")


def __delly_batch(args, file_sets):
    """ Writes batch scripts for running Delly """
    with open("{}/2_SVs/delly_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name delly\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 4\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%12\n".format(len(file_sets.keys())))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.delly.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.delly.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("REF={}\n\n".format(args.r))

        outfile.write('echo "Running delly on ${SAMPLE}"\n\n')

        outfile.write("bash /home.roaming/s4097594/SV_pipeline/scripts/delly.sh -s ${SAMPLE} -r ${REF} -o ${OUTDIR}")


def __freebayes_batch(args, file_sets):
    """ Writes batch scripts for running freebayes """
    with open("{}/2_SVs/freebayes_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name freebayes\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 4\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%12\n".format(len(file_sets.keys())))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.fb.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.fb.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate beatson_py3\n\n")

        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("REF={}\n\n".format(args.r))

        outfile.write('echo "Running freebayes on ${SAMPLE}"\n\n')

        outfile.write("bash /home.roaming/s4097594/SV_pipeline/scripts/freebayes.sh -s ${SAMPLE} -r ${REF} -o ${OUTDIR}")


def __gatk_batch(args, file_sets):
    """ Writes batch scripts for running GATK HaplotypeCaller
    Hard coded GATK threads as was defaulting to 4 threads.
    If encounter an out of memory error, check number of requested threads. """
    with open("{}/2_SVs/gatk_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name gatk\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 8\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%6\n".format(len(file_sets.keys())))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.gatk.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.gatk.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("REF={}\n\n".format(args.r))

        outfile.write('echo "Running gatk on ${SAMPLE}"\n\n')

        outfile.write("bash /home.roaming/s4097594/SV_pipeline/scripts/gatk.sh -s ${SAMPLE} -r ${REF} -o ${OUTDIR}")


def __pindel_batch(args, file_sets):
    """ Writes batch scripts for running Pindel.
    Hard coded Pindel threads as does not use on average more than 6 threads and
    only portions of Pindel code runs on multiple threads. """
    with open("{}/2_SVs/pindel_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name pindel\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 16\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%3\n".format(len(file_sets.keys())))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.pindel.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.pindel.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("REF={}\n".format(args.r))
        outfile.write("THREADS=16\n\n")

        outfile.write('echo "Running pindel on ${SAMPLE}"\n\n')

        outfile.write("bash /home.roaming/s4097594/SV_pipeline/scripts/pindel.sh -t ${THREADS} -s ${SAMPLE} -r ${REF} -o ${OUTDIR}")


def __filter_batch(args, file_sets):
    """ Writes batch script for filtering and intersection data """
    tools = ""
    if args.D:
        tools += "-D "
    if args.F:
        tools += "-F "
    if args.G:
        tools += "-G "
    if args.P:
        tools += "-P "

    with open("{}/2_SVs/filter_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name filter_intersect\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 4\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%12\n".format(len(file_sets.keys())))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.filter_intersect.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.filter_intersect.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")
        outfile.write("OUTDIR={}\n\n".format(args.o))
        outfile.write("REF={}\n".format(args.r))
        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

        outfile.write('echo "Filtering ${SAMPLE}"\n\n')
        outfile.write("bash /home.roaming/s4097594/SV_pipeline/scripts/filter_bcf.sh -s ${{SAMPLE}} -o ${{OUTDIR}} {}".format(tools))

        outfile.write('echo "Getting intersect"\n\n')
        outfile.write("python3 /home.roaming/s4097594/SV_pipeline/scripts/intersect.py -r ${{REF}} -o ${{OUTDIR}} -s ${{SAMPLE}} -S -I {}".format(tools))



def __collate_batch(args, file_sets):
    """ Writes batch script for collating data (NO LONGER USED)"""
    tools = ""
    if args.D:
        tools += "-D "
    if args.F:
        tools += "-F "
    if args.G:
        tools += "-G "
    if args.P:
        tools += "-P "

    with open("{}/2_SVs/collate_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name collate\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 16\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --output={}/batch_output/s%A.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write("REF={}\n".format(args.r))
        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("SAMPLES={}/sample_names.txt\n\n".format(args.o))

        outfile.write('echo "Collating data"\n\n')

        outfile.write("python3 /home.roaming/s4097594/SV_pipeline/scripts/collate.py -r ${{REF}} -o ${{OUTDIR}} -s ${{SAMPLES}} -I {}".format(tools))


def __intersect_batch(args, file_sets):
    """ Writes batch script for getting intersect of SVs for each sample (NOW RUN WITH FILTER STAGE) """
    tools = ""
    if args.D:
        tools += "-D "
    if args.F:
        tools += "-F "
    if args.G:
        tools += "-G "
    if args.P:
        tools += "-P "

    with open("{}/2_SVs/intersect_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name intersect\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 1\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=1-{}%48\n".format(len(file_sets.keys())))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write("REF={}\n".format(args.r))
        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))


        outfile.write('echo "Getting intersect"\n\n')

        outfile.write("python3 /home.roaming/s4097594/SV_pipeline/scripts/intersect.py -r ${{REF}} -o ${{OUTDIR}} -s ${{SAMPLE}} -S -I {}".format(tools))


def __matrix_batch(args, file_sets):
    """ Writes batch script for collating data """
    rec = ""
    if args.E:
        rec += "-E {} ".format(args.E)

    tools = ""
    if args.D:
        tools += "-D "
    if args.F:
        tools += "-F "
    if args.G:
        tools += "-G "
    if args.P:
        tools += "-P "

    with open("{}/2_SVs/matrix_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name matrix\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 16\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --output={}/batch_output/s%A.mtx.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A.mtx.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write("REF={}\n".format(args.r))
        outfile.write("OUTDIR={}\n".format(args.o))
        outfile.write("SAMPLES={}/sample_names.txt\n\n".format(args.o))

        outfile.write('echo "Collating data"\n\n')

        outfile.write("python3 /home.roaming/s4097594/SV_pipeline/scripts/matrix_fast.py -r ${{REF}} -o ${{OUTDIR}} -s ${{SAMPLES}} {}-S -I {}".format(rec, tools))

def __snpEff_batch(args, file_sets):
        with open("{}/2_SVs/snpEff_batch.sh".format(args.o), "w") as outfile:
            outfile.write("#!/bin/bash\n")
            outfile.write("#SBATCH --job-name snpEff\n")
            outfile.write("#SBATCH --nodes 1\n")
            outfile.write("#SBATCH --cpus-per-task 1\n")
            outfile.write("#SBATCH --nice\n")
            outfile.write("#SBATCH --array=1-{}%48\n".format(len(file_sets.keys())))
            outfile.write("#SBATCH --output={}/batch_output/s%A.snpEff.out\n".format(args.o))
            outfile.write("#SBATCH --error={}/batch_output/s%A.snpEff.err\n".format(args.o))
            outfile.write("#SBATCH --time=3-00:00\n\n")

            outfile.write("source activate beatson_py3\n\n")

            outfile.write("OUTDIR={}\n".format(args.o))
            outfile.write('SAMPLE=$( sed -n "${{SLURM_ARRAY_TASK_ID}}p" {}/sample_names.txt )\n\n'.format(args.o))

            outfile.write('echo "Annotating vcf"\n\n')

            outfile.write("snpEff -c /QNAP/caitlin/snpEff_db/snpEff.config EC958 {}/2_SVs/Final/${SAMPLE}.vcf > {}/2_SVs/Annotated/${SAMPLE}.annot.vcf".format(args.o, args.o))



def __samplot_batch(args, file_sets):
    """ Writes batch script for plotting deletions and translocations wtih samplot """
    with open("{}/2_SVs/samplot_batch.sh".format(args.o), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name samplot\n")
        outfile.write("#SBATCH --nodes 1\n")
        outfile.write("#SBATCH --cpus-per-task 1\n")
        outfile.write("#SBATCH --nice\n")
        outfile.write("#SBATCH --array=0-{}%16\n".format(len(file_sets.keys())//30))
        outfile.write("#SBATCH --output={}/batch_output/s%A_%a.samplot.out\n".format(args.o))
        outfile.write("#SBATCH --error={}/batch_output/s%A_%a.samplot.err\n".format(args.o))
        outfile.write("#SBATCH --time=3-00:00\n\n")

        outfile.write("source activate caitlin\n\n")

        outfile.write('echo "Drawing samplots group ${SLURM_ARRAY_TASK_ID}"\n\n')

        outfile.write("python3 /home.roaming/s4097594/SV_pipeline/scripts/draw_samplots.py -a ${{SLURM_ARRAY_TASK_ID}} -o {} -v {}/2_SVs/coregenome.vcf".format(args.o, args.o))
        outfile.write("python3 /home.roaming/s4097594/SV_pipeline/scripts/draw_samplots.py -a ${{SLURM_ARRAY_TASK_ID}} -o {} -v {}/2_SVs/pangenome.vcf".format(args.o, args.o))


def __raxml_batch(args, file_sets):
    """ Writes batch script for running RaxML """
        with open("{}/3_Trees/run_raxml.sh".format(args.o), "w") as outfile:
            outfile.write("#!/bin/bash\n")
            outfile.write("#SBATCH --job-name raxml\n")
            outfile.write("#SBATCH --nodes 1\n")
            outfile.write("#SBATCH --cpus-per-task {}\n".format(args.t))
            outfile.write("#SBATCH --nice\n")
            outfile.write("#SBATCH --output={}/batch_output/s%A_%a.raxml.out\n".format(args.o))
            outfile.write("#SBATCH --error={}/batch_output/s%A_%a.raxml.err\n".format(args.o))
            outfile.write("#SBATCH --time=3-00:00\n\n")

            outfile.write("source activate caitlin\n\n")

            outfile.write("raxmlHPC -f a -N 1000 -s core_snp.phy -x 12345 -p 12345 -m GTRCAT -T {} -n core_snp.nwk".format(args.t))
            outfile.write("raxmlHPC -f a -N 1000 -s core_indel.phy -x 12345 -p 12345 -m BINCAT -T {} -n core_indel.nwk".format(args.t))
            outfile.write("raxmlHPC -f a -N 1000 -s core.phy -q core.partition -x 12345 -p 12345 -m GTRCAT -T {} -n core.nwk".format(args.t))

def __mrbayes_batch(args, file_sets):
    """ Writes batch script for running MrBayes """
    pass

def main():
    args = parse_args()
    file_sets = get_samples(args.s) # get sample names and read file names
    write_batch_files(args, file_sets)

if __name__ == '__main__':
    main()
