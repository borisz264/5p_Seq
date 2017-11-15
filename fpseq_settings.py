import os
import ConfigParser
import simplejson
import shutil
import datetime
import sys

import fpseq_utils

class fpseq_settings:
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings(settings_file)

    def get_settings_file(self):
        return self.settings_file

    def get_property(self, property, default=None):
        try:
            if not property in self.settings and default != None:
                return default
            return self.settings[property]
        except:
            print self.settings
            raise  ValueError('cannot find %s' % property)

    def get_rdir(self):
        fpseq_utils.make_dir(self.rdir)
        return self.rdir

    def iter_lib_settings(self):
        for i in range(len(self.sample_names)):
            yield fpseq_lib_settings(self, self.sample_names[i], self.read1_handles[i], self.read2_handles[i])

    def process_settings(self, settings_file):
        """
        - reads the settings file and converts str to float, list, etc.
        - stores result in self.settings as a dict()
        - CRITICAL NOTE: All keys must be lower case
        """
        int_keys = ['rand_barcode_length', 'ncrna_genomesasparsed', 'ncrna_genomesaindexnbases', 'outfiltermultimapnmax',
                    'alignsjdboverhangmin', 'alignsjoverhangmin', 'genomic_genomesasparsed', 'genomic_genomesaindexnbases', 'max_read_length']
        #float_keys = []
        str_keys = ['star_ncrna_dir', 'ncrna_sequence_dir', 'star_genome_dir', 'genome_sequence_dir', 'annotation_gtf_file']
        #boolean_keys = []
        list_str_keys = ['read1_files', 'read2_files', 'sample_names']
        #list_float_keys = ['concentrations', 'input_rna']
        extant_files = ['genome_sequence_dir', 'annotation_gtf_file']
        config = ConfigParser.ConfigParser()
        config.read(settings_file)
        settings = {}
        for section in config.sections():
            for option in config.options(section):
                settings[option] = config.get(section, option)
                settings[section] = True
        for k in int_keys:
            settings[k] = int(settings[k])
        for k in str_keys:
            settings[k] = settings[k]
        #for k in float_keys:
        #    settings[k] = float(settings[k])
        #for k in boolean_keys:
        #    if not settings[k].lower() in ['true', 'false']:
        #        raise ValueError(
        #          'Boolean value %s must be "true" or "false"' % k)
        #    settings[k] = settings[k].lower() == 'true'
        #for k in list_float_keys:
        #    settings[k] = map(float, simplejson.loads(settings[k]))
        #for k in list_int_keys:
        #    settings[k] = map(int, simplejson.loads(settings[k]))
        for k in list_str_keys:
            settings[k] = simplejson.loads(settings[k])
        self.fqdir = settings['fastq_dir']
        self.sample_names = settings['sample_names']
        #for paired end reads, there are now 2 fastq files per sample
        self.read1_files = settings['read1_files']
        self.read2_files = settings['read2_files']
        self.read1_handles = [os.path.join(self.fqdir, read1_file) if read1_file != '' else None for read1_file in self.read1_files]
        self.read2_handles = [os.path.join(self.fqdir, read2_file) if read2_file != '' else None for read2_file in self.read2_files]
        for file_handle in self.read1_handles+self.read2_handles:
            try:
                assert fpseq_utils.file_exists(file_handle)
            except:
                if file_handle is None:
                    print 'Warning: No input file supplied for a read file, this is only currently supported for Read2'
                else:
                    print 'ERROR: nonexistent file ', file_handle
                    sys.exit()
        for k in extant_files:
            try:
                assert fpseq_utils.file_exists(settings[k])
            except AssertionError:
                print 'file %s does not exist' % settings[k]
                sys.exit()
        self.settings = settings
        self.rdir = settings['results_dir']
        fpseq_utils.make_dir(self.rdir)
        shutil.copy(settings_file, self.rdir)

    def get_log(self):
        log = os.path.join(
          self.get_rdir(),
          'log.txt')
        return log

    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

    ##########################
    # Global File Getters
    ##########################
    def get_star_genome_dir(self):
        index = self.get_property('star_genome_dir')
        return index

    def get_star_ncrna_dir(self):
        index = self.get_property('star_ncrna_dir')
        return index

    def get_genome_sequence_dir(self):
        genome_dir = self.get_property('genome_sequence_dir')
        return genome_dir

    def get_ncrna_sequence_dir(self):
        genome_dir = self.get_property('ncrna_sequence_dir')
        return genome_dir

    def get_genome_sequence_files(self, allowed_endings=['.fa']):
        genome_dir = self.get_property('genome_sequence_dir')
        fasta_files = set()
        for file in os.listdir(genome_dir):
            for ending in allowed_endings:
                if file.endswith(ending):
                    fasta_files.add(os.path.join(self.get_property('genome_sequence_dir'), file))
        return sorted(fasta_files)

    def get_annotation_GTF_file(self):
        anno_file = self.get_property('annotation_gtf_file')
        return anno_file

class fpseq_lib_settings:
    def __init__(self, experiment_settings, sample_name, read1_filehandle, read2_filehandle):
        self.experiment_settings = experiment_settings
        self.sample_name = sample_name
        self.read1_filehandle = read1_filehandle
        self.read2_filehandle = read2_filehandle
    def get_property(self, property):
        return self.experiment_settings.get_property(property)

    def get_log(self):
        fpseq_utils.make_dir(os.path.join(self.experiment_settings.get_rdir(), 'logs'))
        log = os.path.join(
          self.experiment_settings.get_rdir(),
          'logs',
          '%(sample_name)s.log' %
           {'sample_name': self.sample_name})
        return log

    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

    ########################
    #File Getters
    ########################
    def get_read1(self):
        return self.read1_filehandle

    def get_read2(self):
        return self.read2_filehandle

    def get_debarcoded_read1(self):
        debarcoded_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'debarcoded',
          '%(sample_name)s_1.fastq.gz' %
           {'sample_name': self.sample_name})
        return debarcoded_reads

    def get_ncrna_mapping_log(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'ncrna_mapped_reads', '%(sample_name)sLog.final.out' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_ncrna_mapped_reads_prefix(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'ncrna_mapped_reads',
                                    '%(sample_name)s' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_genome_mapped_reads_prefix(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'genome_mapped_reads', '%(sample_name)s' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_ncrna_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'ncrna_mapped_reads', '%(sample_name)sAligned.sortedByCoord.out.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_ncrna_unmapped_read1(self):
        unmapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'ncrna_mapped_reads', '%(sample_name)sUnmapped.out.mate1' % {'sample_name': self.sample_name})
        return unmapped_reads

    def get_ncrna_unmapped_read2(self):
        unmapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'ncrna_mapped_reads', '%(sample_name)sUnmapped.out.mate2' % {'sample_name': self.sample_name})
        return unmapped_reads

    def get_ncrna_most_common_reads(self):
        unmapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'QC', '%(sample_name)s.common_ncrna_fragments.tsv' % {'sample_name': self.sample_name})
        return unmapped_reads

    def get_genome_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'genome_mapped_reads', '%(sample_name)sAligned.sortedByCoord.out.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_transcript_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'genome_mapped_reads', '%(sample_name)sAligned.toTranscriptome.out.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_transcript_counts(self):
        sequence_counts = os.path.join(
          self.experiment_settings.get_rdir(),
          'transcript_counts',
          '%(sample_name)s.counts.pkl' %
           {'sample_name': self.sample_name})
        return sequence_counts


    #####################
    # Checks for existence
    #####################

    def debarcoded_reads_exist(self):
        trimmed_reads = self.get_debarcoded_read1()
        return fpseq_utils.file_exists(trimmed_reads)

    def ncrna_mapped_reads_exist(self):
        mapped_reads = self.get_ncrna_mapped_reads()
        return fpseq_utils.file_exists(mapped_reads)

    def genome_mapped_reads_exist(self):
        mapped_reads = self.get_genome_mapped_reads()
        return fpseq_utils.file_exists(mapped_reads)