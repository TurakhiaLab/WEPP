# rule subsample:
#     input:
#         "build/Makefile",
#     output:
#         "build/wbe"
#     shell:
#         "cd build && make -j"