__author__ = 'boris zinshteyn'
"""
Intended for processing of ribosome footprint profiling data from mammalian cells
"""
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import argparse
import subprocess

import fpseq_settings
import fpseq_utils

class experiment:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.settings.write_to_log('Initializing experiment %s' % self.settings.get_property('experiment_name'))
        self.num_libs = len([x for x in settings.iter_lib_settings()])
        self.move_molecular_barcode_to_name()
        self.make_ncRNA_mapping_index()
        self.make_genome_mapping_index()
        self.map_reads_to_ncrna()
        self.map_reads_to_genome()
        # self.settings.write_to_log('loading genome sequence')
        # self.genome = fpseq_utils.genome_sequence(self.settings.get_genome_sequence_files())
        # self.settings.write_to_log('loading GTF annotations')
        # self.GTF_annotations = fpseq_utils.gtf_data(self.settings.get_annotation_GTF_file())
        # self.initialize_libs()
        self.settings.write_to_log('Finished initializing experiment %s\n' % self.settings.get_property('experiment_name'))

    def move_molecular_barcode_to_name(self):
        if self.settings.get_property('rand_barcode_length') == 0:
            self.settings.write_to_log('rand_barcode_length set to zero, not trimming reads')
            return
        else:
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.debarcoded_reads_exist():
                    break
            else: #else clause executes if break did not occur
                self.settings.write_to_log('using existing trimmed reads')
                return
        self.settings.write_to_log('debarcoding reads')
        fpseq_utils.make_dir(self.rdir_path('debarcoded'))
        fpseq_utils.parmap(lambda lib_setting: self.move_molecular_barcode_to_name_one_lib(lib_setting),
                           self.settings.iter_lib_settings(), nprocs=self.threads)
        self.settings.write_to_log('done debarcoding reads')

    def move_molecular_barcode_to_name_one_lib(self, lib_setting):
        lib_setting.write_to_log('debarcoding reads')
        in_fastq_file_name = lib_setting.get_read1()
        out_fastq_file_name = lib_setting.get_debarcoded_read1()
        num_chars = lib_setting.get_property('rand_barcode_length')
        out_fastq_file = fpseq_utils.aopen(out_fastq_file_name, mode='w')
        for quartet in fpseq_utils.iter4Lines(in_fastq_file_name):
            name_line, seq_line, garbage, q_scores = quartet
            barcode = seq_line[:num_chars]
            out_fastq_file.write('%s %s\n' % (name_line.rstrip('\n'), barcode))
            out_fastq_file.write(seq_line[num_chars:])
            out_fastq_file.write(garbage)
            out_fastq_file.write(q_scores[num_chars:])
        out_fastq_file.close()
        lib_setting.write_to_log('done debarcoding reads')

    def make_ncRNA_mapping_index(self):
        make_index = False
        if fpseq_utils.file_exists(self.settings.get_property('star_ncrna_dir')):
            self.settings.write_to_log('STAR index exists at %s' % self.settings.get_property('star_ncrna_dir'))
            self.settings.write_to_log('using existing STAR index')
        else:
            make_index = True
            fpseq_utils.make_dir(self.settings.get_property('star_ncrna_dir'))
        if make_index:
            self.settings.write_to_log('building STAR index')
            command_to_run = 'STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s*.fa --genomeSAsparseD %d  --genomeSAindexNbases %d 1>>%s 2>>%s' % (
                self.threads, self.settings.get_property('star_ncrna_dir'), self.settings.get_ncrna_sequence_dir(),
                self.settings.get_property('ncrna_genomesasparsed'), self.settings.get_property('ncrna_genomesaindexnbases'),
                self.settings.get_log(), self.settings.get_log())
            self.settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell=True).wait()
        self.settings.write_to_log('STAR index ready')

    def make_genome_mapping_index(self):
        make_index = False
        if fpseq_utils.file_exists(self.settings.get_property('star_genome_dir')):
            self.settings.write_to_log('STAR index exists at %s' % self.settings.get_property('star_genome_dir'))
            self.settings.write_to_log('using existing STAR index')
        else:
            make_index = True
            fpseq_utils.make_dir(self.settings.get_property('star_genome_dir'))
        if make_index:
            self.settings.write_to_log('building STAR index')
            command_to_run = 'STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s*.fa --genomeSAsparseD %d --genomeSAindexNbases %d 1>>%s 2>>%s' % (
                self.threads, self.settings.get_property('star_genome_dir'), self.settings.get_genome_sequence_dir(),
                self.settings.get_property('genomic_genomesasparsed'), self.settings.get_property('genomic_genomesaindexnbases'),
                self.settings.get_log(), self.settings.get_log())
            self.settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell=True).wait()
        self.settings.write_to_log('STAR index ready')

    def remove_adaptor(self):
        if self.settings.get_property('adaptor_3p_sequence') == '':
            self.settings.write_to_log('adaptor_3p_sequence not provided, not trimming adaptors fromreads')
            return
        elif not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                self.settings.write_to_log('using existing adaptorless reads')
                return
        else:
            self.settings.write_to_log('adaptor removal forced')
        self.settings.write_to_log('removing adaptors')
        fpseq_utils.make_dir(self.rdir_path('adaptor_removed'))
        map(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting, self.threads),
                          self.settings.iter_lib_settings())
        self.settings.write_to_log('done removing adaptors')

    def remove_adaptor_one_lib(self, lib_settings, threads):
        lib_settings.write_to_log('adaptor trimming')
        """
        -x specifies the 3' adaptor to trim from the forward read
        -Q specifies the lowest acceptable mean read quality before trimming
        -l specifies the minimum post-trimming read length
        -L specifies the maximum post-trimming read length
        -o is the output prefix
        --threads specifies number of threads to use
        """

        if self.settings.get_property('trim_5p') == 0:
            reads_to_use = lib_settings.get_fastq_gz_file()
        else:
            reads_to_use = lib_settings.get_debarcoded_read1()

        command_to_run = 'skewer -x %s -Q %d -l %d -L %d -o %s --quiet --threads %s %s 1>>%s 2>>%s' % (
            self.settings.get_property('adaptor_3p_sequence'),
            self.settings.get_property('sequence_quality_cutoff'),
            self.settings.get_property('min_post_trimming_length'),
            self.settings.get_property('max_post_trimming_length'),
            lib_settings.get_adaptor_trimmed_reads(prefix_only=True),
            threads,
            reads_to_use,
            lib_settings.get_log(), lib_settings.get_log())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        compression_command = 'gzip %s-trimmed.fastq' % (lib_settings.get_adaptor_trimmed_reads(prefix_only=True))
        lib_settings.write_to_log(compression_command)
        subprocess.Popen(compression_command, shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

    def map_reads_to_ncrna(self):
        """
        map all reads using STAR and discard any that map to rRNA
        :return:
        """
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.ncrna_mapped_reads_exist():
                break
        else:
            self.settings.write_to_log('using existing noncoding RNA mapped reads')
            return
        self.settings.write_to_log('mapping reads')
        fpseq_utils.make_dir(self.rdir_path('ncrna_mapped_reads'))
        map(lambda lib_setting: self.map_one_library_to_ncrna(lib_setting, self.threads), self.settings.iter_lib_settings())
        self.settings.write_to_log( 'finished mapping reads to noncoding RNA')

    def map_one_library_to_ncrna(self, lib_settings, threads):
        lib_settings.write_to_log('mapping reads to ncRNA')
        if lib_settings.get_read2() is not None :
            command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 8000000000 --genomeDir %s --readFilesIn %s %s --readFilesCommand gunzip -c ' \
                             '--outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax %d --outWigType wiggle read1_5p --outFileNamePrefix %s ' \
                             '--outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                             (threads, self.settings.get_star_ncrna_dir(),
                              lib_settings.get_debarcoded_read1(), lib_settings.get_read2(),
                              self.settings.get_property('outfiltermultimapnmax'),
                              lib_settings.get_ncrna_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        else:
            command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 8000000000 --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c ' \
                             '--outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax %d --outWigType wiggle read1_5p --outFileNamePrefix %s ' \
                             '--outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                             (threads, self.settings.get_star_ncrna_dir(),
                              lib_settings.get_debarcoded_read1(),
                              self.settings.get_property('outfiltermultimapnmax'),
                              lib_settings.get_ncrna_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        command_to_run = 'samtools index %s' % (lib_settings.get_ncrna_mapped_reads())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('done mapping reads to ncRNA')

    def map_reads_to_genome(self):
        """
        map all reads using STAR
        :return:
        """
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.genome_mapped_reads_exist():
                break
        else:
            self.settings.write_to_log('using existing genome-mapped reads')
            return
        self.settings.write_to_log('mapping reads to genome')
        fpseq_utils.make_dir(self.rdir_path('genome_mapped_reads'))
        map(lambda lib_setting: self.map_one_library_to_genome(lib_setting, self.threads), self.settings.iter_lib_settings())
        self.settings.write_to_log( 'finished mapping reads to genome')

    def map_one_library_to_genome(self, lib_settings, threads):
        lib_settings.write_to_log('mapping_reads')
        if lib_settings.get_read2() is not None :
            command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 8000000000 --genomeDir %s --readFilesIn %s %s ' \
                             '--outSAMtype BAM SortedByCoordinate --sjdbGTFfile %s --sjdbOverhang %d --alignSJDBoverhangMin %d --alignSJoverhangMin %d ' \
                             '--outFilterType BySJout --outFilterMultimapNmax %d --outWigType wiggle read1_5p --outFileNamePrefix %s' \
                             ' --outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                             (threads, self.settings.get_star_genome_dir(),
                              lib_settings.get_ncrna_unmapped_read1(), lib_settings.get_ncrna_unmapped_read2(),
                              self.settings.get_annotation_GTF_file(),
                              self.settings.get_property('max_read_length')-1,
                              self.settings.get_property('alignsjdboverhangmin'),
                              self.settings.get_property('alignsjoverhangmin'),
                              self.settings.get_property('outfiltermultimapnmax'),
                              lib_settings.get_genome_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        else:
            command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 8000000000 --genomeDir %s --readFilesIn %s ' \
                             '--outSAMtype BAM SortedByCoordinate --sjdbGTFfile %s --sjdbOverhang %d --alignSJDBoverhangMin %d --alignSJoverhangMin %d ' \
                             '--outFilterType BySJout --outFilterMultimapNmax %d --outWigType wiggle read1_5p --outFileNamePrefix %s' \
                             ' --outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                             (threads, self.settings.get_star_genome_dir(),
                              lib_settings.get_ncrna_unmapped_read1(),
                              self.settings.get_annotation_GTF_file(),
                              self.settings.get_property('max_read_length')-1,
                              self.settings.get_property('alignsjdboverhangmin'),
                              self.settings.get_property('alignsjoverhangmin'),
                              self.settings.get_property('outfiltermultimapnmax'),
                              lib_settings.get_genome_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        command_to_run = 'samtools index %s' % (lib_settings.get_genome_mapped_reads())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('mapping_reads done')

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        fpseq_utils.make_dir(self.rdir_path('transcript_counts'))
        self.libs = []
        fpseq_utils.parmap(lambda lib_settings: ribo_lib.assign_tx_reads(self, self.settings, lib_settings), self.settings.iter_lib_settings(), nprocs = self.threads)
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')

    def find_lib_by_sample_name(self, sample_name):
        for lib in self.libs:
            if lib.lib_settings.sample_name == sample_name:
                return lib
        assert False #if this triggers, your settings file is broken.

    def initialize_lib(self, lib_settings):
        lib = ribo_lib.ribo_lib(self.settings, lib_settings)
        self.libs.append(lib)

    def make_tables(self):
        fpseq_utils.make_dir(self.rdir_path('tables'))
        ribo_tables.make_readthrough_table(self)
        ribo_tables.make_detailed_readthrough_table(self)
        #ribo_tables.transcriptome_features_table(self)
        ribo_tables.make_cds_rpkm_table(self)
        ribo_tables.make_cds_counts_table(self)

    def make_plots(self):

        fpseq_utils.make_dir(self.rdir_path('plots'))
        ribo_plotting.plot_fragment_length_distributions(self)
        ribo_plotting.plot_frame_distributions(self)
        ribo_plotting.plot_start_codon_average(self, up = 60, down = 180)
        ribo_plotting.plot_stop_codon_average(self, up = 180, down = 60)
        ribo_plotting.plot_stop_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')
        ribo_plotting.plot_start_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')
        ribo_plotting.plot_readthrough_box(self)
        ribo_plotting.plot_readthrough_box(self, log=True)

    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def make_counts_table(self, fractional=False):
        """
        write out number of fragments mapping to each TL in each dataset
        :param fractional: if True, replace raw counts with library fraction in reads per million
        :return:
        """
        if fractional:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'rpm.txt'), 'w')
        else:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'raw_counts.txt'), 'w')

        header = 'sequence name\t' + '\t'.join([lib.lib_settings.sample_name for lib in self.libs]) + '\n'
        summary_file.write(header)
        if fractional:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' % ((10**6)*lib.pool_sequence_mappings[sequence_name].fragment_count/float(lib.total_mapped_fragments)) for lib in self.libs]))
                summary_file.write(out_line)
        else:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' %
                                                    lib.pool_sequence_mappings[sequence_name].fragment_count
                                                    for lib in self.libs]))
                summary_file.write(out_line)
        summary_file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--make-tables",
                        help="Makes tables.",
                        action='store_true')
    parser.add_argument("--perform-qc",
                        help="performs quality control analysis.",
                        action='store_true')
    parser.add_argument("--make-plots",
                        help="Makes plots.",
                        action='store_true')
    parser.add_argument("--all-tasks",
                        help="Makes plots, tables",
                        action='store_true')
    parser.add_argument("--threads",
                        help="Max number of processes to use",
                        type = int, default = 8)
    args = parser.parse_args()

    return args

def main():
    """
    """
    args = parse_args()
    settings = fpseq_settings.fpseq_settings(args.settings_file)
    ribo_experiment = experiment(settings, args.threads)
    print 'experiment ready'

    if args.perform_qc or args.all_tasks:
        settings.write_to_log('performing QC')
        ribo_experiment.perform_qc()
        settings.write_to_log('done performing QC')

    if args.make_tables or args.all_tasks:
        print 'tables'
        settings.write_to_log('making tables')
        ribo_experiment.make_tables()
        settings.write_to_log('done making tables')

    if args.make_plots or args.all_tasks:
        print 'plots'
        settings.write_to_log('making plots')
        ribo_experiment.make_plots()
        settings.write_to_log('done making plots')

if __name__ == "__main__":
    # execute only if run as a script
    main()