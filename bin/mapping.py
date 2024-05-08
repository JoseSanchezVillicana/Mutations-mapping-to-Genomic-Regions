import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gc

def load_gff(file_path: str) -> pd.DataFrame:
    columns = ["seqnames", "feature", "starts", "ends"]
    return pd.read_csv(
        filepath_or_buffer= file_path,
        names= columns,
        header= None,
        sep= '\t',
        comment= '#'
    )

def load_vcf(file_path: str) -> pd.DataFrame:
    columns = ["seqnames", "pos", "id", 'REF', 'ALT']
    return pd.read_csv(
        filepath_or_buffer= file_path,
        names= columns,
        header= None,
        sep= '\t',
        comment= '#',
        dtype= {
            'seqnames': 'string',
            'pos': 'int64',
            'id': 'int64',
            'REF': 'string',
            'ALT' : 'string'
        }
    )

def binary_search(pos: int, genome_regions: pd.DataFrame) -> str:
    search_pool = pd.Index(genome_regions)
    i = 0
    f = len(search_pool) - 1
    while i <= f:
        m = (i + f) // 2
        if pos >= search_pool[m][2] and pos <= search_pool[m][3]:
            return search_pool[m][1]
        elif pos > search_pool[m][3]: #Means its in right region of genome
            i = m + 1
        else:
            f = m - 1
    return 'unmapped'

def get_intergenic_regions(genes_df: pd.DataFrame) -> pd.DataFrame:
    individual_intergenic_regions = []
    for i in range(1, len(genes_df)):
        if genes_df.iloc[i]["seqnames"] == genes_df.iloc[i-1]["seqnames"]:
            intergenic = {
                "seqnames": genes_df.iloc[i - 1]["seqnames"],
                "starts": genes_df.iloc[i - 1]["ends"] + 1,
                "ends": genes_df.iloc[i]["starts"] - 1,
                "feature": "intergenic"
            }
            if intergenic['ends'] >= intergenic['starts']:
                individual_intergenic_regions.append(intergenic)

    intergenic_regions = {
        'seqnames' : [],
        'feature' : [],
        'starts' : [],
        'ends' : []
    }
    for intergen_region in individual_intergenic_regions:
        for feature, value in intergen_region.items():
            intergenic_regions[feature].append(value)
    
    return pd.DataFrame(intergenic_regions)

def get_intron_regions(coding_df: pd.DataFrame) -> pd.DataFrame:
    individual_intronic_regions = []

    for i in range(len(coding_df) - 1):
        if coding_df.iloc[i]["feature"] == 'gene' or coding_df.iloc[i]["feature"] == 'exon' or coding_df.iloc[i]["feature"] == 'CDS':
            if coding_df.iloc[i + 1]["feature"] == 'exon' or coding_df.iloc[i + 1]["feature"] == 'CDS' or coding_df.iloc[i + 1]["feature"] == 'intergenic':
                if coding_df.iloc[i + 1]['starts'] > coding_df.iloc[i]['ends']:
                    intronic = {
                        "seqnames": coding_df.iloc[i]['seqnames'],
                        "starts": coding_df.iloc[i]["ends"] + 1,
                        "ends": coding_df.iloc[i + 1]["starts"] - 1,
                        "feature": "intron"
                    }
                    individual_intronic_regions.append(intronic)

    intronic_regions = {
        'seqnames' : [],
        'feature' : [],
        'starts' : [],
        'ends' : []
    }
    for intron_region in individual_intronic_regions:
        for feature, value in intron_region.items():
            intronic_regions[feature].append(value)
    
    return pd.DataFrame(intronic_regions)


def get_genomic_regions(reference_genome: str, mutations: str, verbose= None):
    refseqs = {
        '1': 'NC_000001.11',
        '2': 'NC_000002.12',
        '3' : 'NC_000003.12',
        '4' : 'NC_000004.12',
        '5' : 'NC_000005.10',
        '6' : 'NC_000006.12',
        '7' : 'NC_000007.17',
        '8' : 'NC_000008.11',
        '9' : 'NC_000009.12',
        '10' : 'NC_000010.11',
        '11' : 'NC_000011.10',
        '12' : 'NC_000012.12',
        '13' : 'NC_000013.11',
        '14' :'NC_000014.9',
        '15' : 'NC_000015.10',
        '16' : 'NC_000016.10',
        '17' : 'NC_000017.11',
        '18' : 'NC_000018.10',
        '19' : 'NC_000019.10',
        '20' : 'NC_000020.11',
        '21' : 'NC_000021.9',
        '22' : 'NC_000022.11',
        'X' : 'NC_000023.11',
        'Y' : 'NC_000024.10',
        'MT' : 'NC_012920.1'
    }

    if verbose: 
        print('\tLoading reference genome...\n')
    gff_df = load_gff(reference_genome)
    mutations_df = load_vcf(mutations).sort_values(by=['pos']).reset_index(drop= True)

    fig, axs = plt.subplots(5, 5, figsize=(15, 15))

    for i, chr in enumerate(refseqs.keys()):
        if chr == '21':
            if verbose:
                print(f'\t----> Processing chromosome {chr}...\n')
            
            chr_gff = gff_df[gff_df["seqnames"] == refseqs[chr]]
            chr_mutations = mutations_df[mutations_df['seqnames'] == chr]
            mutations_positions = list(chr_mutations['pos'])

            # Filter only genes regions
            genes_df = chr_gff[chr_gff["feature"] == "gene"]
            genes_df = genes_df.sort_values(by=["starts", 'ends']).reset_index(drop= True)

            # Filter coding regions
            coding_df = chr_gff[(chr_gff["feature"] == "gene") | (chr_gff['feature'] == 'exon') | (chr_gff['feature'] == 'CDS')]

            # Get intergenic regions
            if verbose:
                print('\tCalculating intergenic regions...\n')
            intergenic_regions = get_intergenic_regions(genes_df= genes_df)
            temp_concatenated = pd.concat([coding_df, intergenic_regions], ignore_index=True).sort_values(by=["starts"]).reset_index(drop= True)

            # Get intronic regions
            if verbose:
                print('\tCalculating intronic regions...\n')
            intronic_regions = get_intron_regions(coding_df= temp_concatenated)
            concatenated_df = pd.concat([intronic_regions, temp_conca], ignore_index=True).sort_values(by=["starts", 'ends']).reset_index(drop= True)
            
            # Mapping regions
            counts = {
                'intergenic' : 0,
                'gene' : 0,
                'unmapped' : 0,
                'CDS' : 0,
                'exon' : 0,
                'intron' : 0
            }
            if verbose:
                print('\tMapping mutations...\n')
            for j, pos in enumerate(mutations_positions):
                if j <= 500:
                    if verbose:
                        print(f'Mapping mutation {j + 1}/{len(mutations_positions)} ...')
                    counts[binary_search(pos, concatenated_df)] += 1

            # Deletion of unnecessary objects
            if verbose:
                print('\tCollecting garbage...\n')
            del chr_gff
            del chr_mutations
            del genes_df
            del intergenic_regions
            del concatenated_df
            gc.collect()

            if verbose:
                print('\tPlotting counts...\n')
            # Plotting
            x = i // 5
            y = i  % 5
            plt.sca(axs[x, y])

            labels = counts.keys()
            sizes = counts.values()
            total = sum(counts.values())
            non_zero_sizes = [size for size in sizes if ((size * 100) / total) > 3]
            non_zero_labels = [label for label, size in zip(labels, sizes) if ((size * 100) / total) > 3]
            colors = cm.tab10.colors 
            plt.pie(
                non_zero_sizes,
                labels= non_zero_labels,
                colors= colors,
                autopct= '%1.2f%%',
                startangle= 90
            )
            plt.title(f'Chromosome {chr}')
            plt.tight_layout()

    return fig

if __name__ == "__main__":
    archivo_gtf = "./../data/GRCH38.p14_trimmed.gff"
    archivo_vcf = './../data/clinvar_trimmed.vcf'
    get_genomic_regions(
        reference_genome= archivo_gtf,
        mutations= archivo_vcf,
        verbose= True
    ).savefig(
        fname= './../results/prueba.png',
        transparent=True,
        bbox_inches='tight'
    );
