import roh_plotting
import variant_plotting
import coverage_plotting
import copy_number_plotting

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt


'''
plot a chromosome on axes
'''
def plot_chrom_on_axes(copy_number, roh, germline_calls, somatic_calls,
                       tumour_coverage, normal_coverage, chrom, axes):

    prepped_copy_number = copy_number_plotting.prepare_at_chrom(copy_number, chrom)
    prepped_roh = roh_plotting.prepare_at_chrom(roh, chrom)
    prepped_somatic_calls = variant_plotting.prepare_at_chrom(somatic_calls, chrom)
    prepped_germline_calls = variant_plotting.prepare_at_chrom(germline_calls, chrom)
    prepped_tumour_coverage = coverage_plotting.prepare_at_chrom(tumour_coverage, chrom, bin=False)
    prepped_normal_coverage = coverage_plotting.prepare_at_chrom(normal_coverage, chrom, bin=True)


    coverage_ylim_MAX = prepped_tumour_coverage.coverage.max() + 10
    coverage_ylim_MIN = prepped_normal_coverage.coverage.min() - 10

    anno_genes = copy_number_plotting.get_gene_annotation_data(chrom)

    axes[0] = variant_plotting.plot(prepped_somatic_calls, axes[0], "somatic")
    axes[1] = variant_plotting.plot(prepped_germline_calls, axes[1], "germline")
    axes[2] = coverage_plotting.plot(prepped_tumour_coverage,
                                     coverage_ylim_MIN, coverage_ylim_MAX, axes[2], "tumour")
    axes[3] = coverage_plotting.plot(prepped_normal_coverage,
                                     coverage_ylim_MIN, coverage_ylim_MAX, axes[3], "normal")
    axes[4] = roh_plotting.plot(prepped_roh, axes[4])
    axes[5] = copy_number_plotting.plot(prepped_copy_number, anno_genes, axes[5])

    return axes

'''
rasterize aces
'''
def rasturize_axes(axes):
    for ax in axes:
        ax.set_rasterized(True)
    return axes


'''
get rid of x axes, labels
force all plots to use same x axis 
'''
def format_axes(ref_axes, axes):
    for ax in axes:
        ax.set_xlim(*ref_axes)
        ax.get_xaxis().set_visible(False)

    return axes


'''
make the plot
'''
def QC_plot(copy_number, roh, germline_calls, somatic_calls,
         tumour_coverage, normal_coverage, pdf):

    #read things in
    copy_number = copy_number_plotting.read(copy_number)
    roh = roh_plotting.read(roh)
    somatic_calls = variant_plotting.read(somatic_calls)
    germline_calls = variant_plotting.read(germline_calls)
    tumour_coverage = coverage_plotting.read(tumour_coverage)
    normal_coverage = coverage_plotting.read(normal_coverage)

    #get chroms to plot on and make sure all data contains those chroms
    chroms = copy_number.Chr.unique()

    for chrom in chroms:

        #make a figure
        fig, axes = plt.subplots(6, gridspec_kw={'height_ratios': [1, 1, 1, 1, 1, 3]},
                                 figsize=(15, 10))

        axes = plot_chrom_on_axes(copy_number, roh, germline_calls,
                                  somatic_calls, tumour_coverage,
                                  normal_coverage, chrom, axes)

        ref_axes = axes[5].get_xlim()
        axes = format_axes(ref_axes, axes)
        axes[5].get_xaxis().set_visible(True)

        axes[0].set_title("({}) Chromosome: {}".format("GRCh37", str(chrom)))
        axes[5].set_xlabel("Position (Mb)")

        axes = rasturize_axes(axes)

        #write out
        pdf.savefig(fig)


def main():

    #testdata
    cn_7 = "/Users/abramsd/work/DATA/QC/titan/Sample007_titan_markers.csv.gz"
    roh_7 = "/Users/abramsd/work/DATA/QC/roh/007/007_roh_ST.txt"
    germ_7 = "/Users/abramsd/work/DATA/QC/germline/007/Sample_007_samtools_germline.vcf"
    som_7 = "/Users/abramsd/work/DATA/QC/somatic/som007.vcf"
    cov_7 = "/Users/abramsd/work/DATA/QC/coverage/merged007_TT"

    cov_norm_7 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/7.bins"


    pdf = matplotlib.backends.backend_pdf.PdfPages("TEST.pdf")
    QC_plot(cn_7, roh_7, germ_7, som_7, cov_7, cov_norm_7, pdf)
    pdf.close()


if __name__ == "__main__":
    main()
