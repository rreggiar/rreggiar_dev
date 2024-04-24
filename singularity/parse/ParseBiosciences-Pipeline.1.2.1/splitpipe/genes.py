#
#   Copyright (c) 2024 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#
# Genome / gene stuff

import pandas as pd
from collections import defaultdict

LOCAL_IMPORT = 0

if LOCAL_IMPORT:
    import utils
else:
    from splitpipe import utils


### ----------------------------------------- Genome stuff -----------------------------------------
def get_genome_specs(gene_info):
    """Get specifications (not data) for genome (e.g. for json)

    Return dict
    """
    new_dict = {}
    all_genes = gene_info["all_genes"]
    genomes = gene_info["genome_list"]
    gen_sizes = []
    for geno in genomes:
        num = len(all_genes[all_genes["genome"] == geno])
        gen_sizes.append(str(num))

    new_dict["genome_dir"] = gene_info.get("genome_dir")
    new_dict["genome_num"] = len(genomes)
    new_dict["genome_names"] = ",".join(genomes)
    new_dict["genome_list"] = genomes
    new_dict["genome_sizes"] = gen_sizes
    new_dict["genome_total_size"] = len(all_genes)

    return new_dict


def report_genome_specs_str(spipe, pref="# Genome_specs "):
    """Get lines reporting genome specifications

    Return str
    """
    pdict, source = spipe.get_runproc_info(mkref=True)
    if not pdict:
        return ""

    str_rows = []
    run_version = pdict["header"]["version"]
    run_date = pdict["header"]["date"]
    str_rows.append(f"{pref}mkref version {run_version}")
    str_rows.append(f"{pref}mkref date   {run_date}")

    geno_dir = pdict["genome_dir"]
    str_rows.append(f"{pref}genome_dir   {geno_dir}")

    genome_list = pdict["genome"]["genome_list"]
    genome_sizes = pdict["genome"]["genome_sizes"]
    in_dict_list = pdict["mkref"]["input_files"]
    for genome, size, in_dict in zip(genome_list, genome_sizes, in_dict_list):
        str_rows.append(f"{pref}number_genes {genome}  {size}")
        str_rows.append(
            f"{pref}in_fasta     {genome}  {in_dict['fasta'].split(',')[0]}"
        )
        str_rows.append(f"{pref}in_gtf       {genome}  {in_dict['gtf'].split(',')[0]}")
    if len(genome_list) > 1:
        str_rows.append(f"{pref}total_genes  {pdict['genome']['genome_total_size']}")

    return "\n".join(str_rows)


def write_gene_info(gene_info, fpath, verb=True, indent=None):
    """Write out contents of gene_info structure to json file

    Return filename
    """
    # Some items need to be "simplifed" for json
    # Make copy of top dict then replace as needed
    new_dict = {}
    for k, v in gene_info.items():
        new_dict[k] = v
    gene_info = new_dict
    # Convert contained int to actual ints (key and val)
    gene_info["genes_to_exons"] = clean_dic_sub_int_dicts(gene_info["genes_to_exons"])
    # Convert dataframe to json str
    if "all_genes" in gene_info:
        gene_info["all_genes"] = gene_info["all_genes"].to_json()

    # Save
    fpath = utils.write_json(gene_info, fpath, indent=indent)
    if verb:
        print(f"Wrote {fpath}")
    return fpath


def clean_dic_sub_int_dicts(d_dic):
    """Replace all dicts int dict with new int key,val dicts

    Return dict
    """
    new_dict = {}
    for k, v in d_dic.items():
        c_dict = clean_int_dict(v)
        new_dict[k] = c_dict
    return new_dict


def clean_int_dict(dic):
    """Create a dict with integet keys and values from given dict

    Return dict
    """
    c_dict = {}
    for k, v in dic.items():
        ki = int(k)
        vi = int(v)
        c_dict[ki] = vi
    return c_dict


def load_gene_info(fpath):
    """Load and prepare gene_info structure from json file

    Return gene_info
    """
    gene_info = utils.read_json(fpath)
    # Should contain simple ints (both key and val)
    gene_info["genes_to_exons"] = clean_dic_sub_int_dicts(gene_info["genes_to_exons"])

    # These should be default dicts
    gene_info["gene_bins"] = defaultdict(list, gene_info["gene_bins"])
    gene_info["genes_to_exons"] = defaultdict(dict, gene_info["genes_to_exons"])

    # Convert json (str) to dataframe if all_genes is here
    if "all_genes" in gene_info:
        gene_info["all_genes"] = pd.io.json.read_json(gene_info["all_genes"])
    return gene_info


# Columns for gtf files / dataframe

GTF_COLS = "Chromosome,Source,Feature,Start,End,Score,Strand,Frame,Attributes".split(
    ","
)


# ---------------------------------------------------------------------------
class GtfData:
    """Class for gtf 'gene' info"""

    def __init__(
        self, data, name="", filename="", gtf_info=None, guess_min=1, standardize=True
    ):
        """Initialize structure

        data = input data; Filename string or dataframe
        """
        self.name = name
        self.filename = filename
        self.flavor = self.att_biotype = self.feat_gene = self.feat_exon = ""
        self.df = None
        # Dataframe or filename?
        if isinstance(data, pd.DataFrame):
            self.df = data
        elif isinstance(data, str):
            # Read as tab delimited strings
            self.df = pd.read_csv(
                data, sep="\t", comment="#", names=GTF_COLS, dtype="object"
            )
            self.filename = data
        else:
            story = f"Given data of unknown type '{type(data)}'"
            raise ValueError(story)

        # Cast coords as int
        self.df["Start"] = self.df["Start"].astype(int)
        self.df["End"] = self.df["End"].astype(int)

        # Need guessed format flavor fields
        if not gtf_info:
            gtf_info, gtf_story = gtf_guess_flavor(data, guess_min)
            if not gtf_info:
                story = "Problem with tgf format flavor guessing...\n" + gtf_story
                raise ValueError(story)

        # Unpack info fields
        self.flavor = gtf_info["flavor"]
        self.at_biotype = gtf_info["at_biotype"]
        self.feat_gene = gtf_info["feat_gene"]
        self.feat_exon = gtf_info["feat_exon"]

        # Standardize format?
        if standardize:
            self._standardize()

    def __del__(self):
        pass

    def __repr__(self):
        ostring = "GtfData (gtf gene info)\n"
        ostring += f"Name:    {self.get_name()}, format flavor {self.get_flavor()}\n"
        ostring += f"Data:    {self.num_records()} records\n"
        return ostring

    def _standardize(self):
        """Standardize key values ..."""
        # Set 'type' to 'biotype'
        if self.at_biotype == "gene_type":
            self.df["Attributes"] = self.df["Attributes"].str.replace(
                "gene_type", "gene_biotype"
            )
            self.att_biotype = "gene_biotype"

    def get_name(self):
        return self.name

    def get_filename(self):
        return self.filename

    def get_flavor(self):
        return self.flavor

    def get_df(self):
        return self.df

    def num_records(self):
        n = 0
        df = self.get_df()
        if df is not None:
            n = len(df)
        return n

    def merge_df(self, data):
        """Merge data into current object"""
        filename = "df"
        # Get dataframe
        if isinstance(data, pd.DataFrame):
            df = data
        elif isinstance(data, GtfData):
            df = data.get_df()
            filename = data.get_filename()
        else:
            story = f"Given data of unknown type '{type(data)}'"
            raise ValueError(story)

        # Merge self dataframe with given one
        cur_df = self.get_df()
        self.df = pd.concat([cur_df, df], axis=0)
        # Concat filename
        self.filename += f",{filename}"


def gtf_guess_flavor(data, guess_min, verb=True):
    """Get specific flavor of gtf in dataframe or dict of these

    data = gtf filename or dataframe with gtf cols
    guess_min = min counts for guessing things

    Return tuple (dict, story)
    """
    story = ""
    if isinstance(data, pd.DataFrame):
        df = data
    elif isinstance(data, str):
        if verb:
            print(f"# Reading gtf from {data}")
        try:
            df = pd.read_csv(
                data, sep="\t", comment="#", names=GTF_COLS, dtype="object"
            )
        except Exception as e:
            print(f"# Problem getting gtf: {e} (exception)")
            return None, story
    else:
        story = f"gtf_guess_flavor given unknown data type: '{type(data)}'"
        raise ValueError(story)

    guess_min = int(guess_min)
    guess = feat_gene = feat_exon = at_biotype = ""
    trust = False

    # Counts to decide on
    if verb:
        print(f"# Tallying data values from {len(df)} records")
    n_feat_gene = df["Feature"].str.count("gene").sum()
    n_feat_exon = df["Feature"].str.count("exon").sum()
    n_feat_cds = df["Feature"].str.count("CDS").sum()
    n_at_gene_biotype = df["Attributes"].str.count("gene_biotype").sum()
    n_at_gene_type = df["Attributes"].str.count("gene_type").sum()

    # Build a string reporting checked value counts
    story = "Tallied the following gtf values\n"
    story = f"  (Min value to count {guess_min})\n"
    story += f"Feature   'gene': {n_feat_gene}\n"
    story += f"Feature   'exon': {n_feat_exon}\n"
    story += f"Feature   'CDS': {n_feat_cds}\n"
    story += f"Attribute 'gene_biotype': {n_at_gene_biotype}\n"
    story += f"Attribute 'gene_type': {n_at_gene_type}\n"

    # Ensembl has 'gene' and 'exon' features and 'gene_biotype' attributes
    if (
        (n_feat_gene >= guess_min)
        and (n_feat_exon >= guess_min)
        and (n_at_gene_biotype >= guess_min)
    ):
        trust = True
        guess = "Ensembl"
        at_biotype = "gene_biotype"
        feat_gene = "gene"
        feat_exon = "exon"
    # Genecode has 'gene' and 'exon' features and 'gene_type' attributes
    elif (
        (n_feat_gene >= guess_min)
        and (n_feat_exon >= guess_min)
        and (n_at_gene_type >= guess_min)
    ):
        trust = True
        guess = "Gencode"
        at_biotype = "gene_type"
        feat_gene = "gene"
        feat_exon = "exon"
    # UCSC has 'CDS' features ... but don't trust this; Keep trust False
    else:
        guess = "???; Not Ensembl or Gencode"

    story += f"Guessed flavor: {guess}\n"

    gtf_info = {
        "trust": trust,
        "flavor": guess,
        "at_biotype": at_biotype,
        "feat_gene": feat_gene,
        "feat_exon": feat_exon,
    }
    return gtf_info, story
