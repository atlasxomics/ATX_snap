import snapatac2 as snap

if __name__ == '__main__':

    fragment_file = snap.datasets.pbmc500(downsample=True)
    data = snap.pp.import_data(
        [fragment_file],
        chrom_sizes=snap.genome.hg38,
        file=['test.h5ad'],
        sorted_by_barcode=False,
    )
    data[0].obs['barcode'] = data[0].obs_names

    print(data)

# fragment_file = "test2.tsv.gz"
# demo = snap.pp.import_data(
#     [fragment_file],
#     chrom_sizes=snap.genome.hg38,
#     file=['test2.h5ad'],
#     sorted_by_barcode=False,
# )
# demo.obs['barcode'] = demo.obs_names

# print(demo)

# # fragment_files = ["test1.tsv.gz", "test2.tsv.gz"]
# # tests = snap.pp.import_data(
# #     fragment_files,
# #     chrom_sizes=snap.genome.hg38,
# #     file=['t1.h5ad', "t2.h5ad"],
# #     sorted_by_barcode=False,
# #     n_jobs=2
# # )
# # print(tests)
# # # for test in tests:
# # #     test.obs["barcode"] = test.obs_names
