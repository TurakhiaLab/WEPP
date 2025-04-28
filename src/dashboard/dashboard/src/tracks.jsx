const tracks = [
  {
    type: "FeatureTrack",
    trackId: "nextstrain-annotations",
    name: "Genes",
    assemblyNames: ['NC_045512'],
    category: [],
    adapter: {
      type: "FromConfigAdapter",
      features: [
        {
          refName: 'NC_045512.2',
          name: "E",
          uniqueId: 4,
          start: 26244,
          end: 26472,
          fill: "#D9AD3D",
        },
        {
          refName: 'NC_045512.2',
          name: "M",
          uniqueId: 5,
          start: 26522,
          end: 27191,
          fill: "#5097BA",
        },
        {
          refName: 'NC_045512.2',
          name: "N",
          uniqueId: 10,
          start: 28273,
          end: 29533,
          fill: "#E67030",
        },
        {
          refName: 'NC_045512.2',
          name: "Orf1a",
          uniqueId: 0,
          start: 265,
          end: 13468,
          fill: "#8EBC66",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF1b",
          uniqueId: 1,
          start: 13467,
          end: 21555,
          fill: "#E59637",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF3a",
          uniqueId: 3,
          start: 25392,
          end: 26220,
          fill: "#AABD52",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF6",
          uniqueId: 6,
          start: 27201,
          end: 27387,
          fill: "#DF4327",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF7a",
          uniqueId: 7,
          start: 27393,
          end: 27759,
          fill: "#C4B945",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF7b",
          uniqueId: 8,
          start: 27755,
          end: 27887,
          fill: "#75B681",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF8",
          uniqueId: 9,
          start: 27893,
          end: 28259,
          fill: "#60AA9E",
        },
        {
          refName: 'NC_045512.2',
          name: "ORF9b",
          uniqueId: 11,
          start: 28283,
          end: 28577,
          fill: "#D9AD3D",
        },
        {
          refName: 'NC_045512.2',
          name: "S",
          uniqueId: 2,
          start: 21562,
          end: 25384,
          fill: "#5097BA",
        },
      ],
    },
    displays: [
      {
        type: "LinearBasicDisplay",
        displayId: "nextstrain-color-display",
        renderer: {
          type: "SvgFeatureRenderer",
          color1: "jexl:get(feature,'fill') || 'black'",
        },
      },
    ],
  },
  // {
  //   type: "DashboardTrack",
  //   trackId: "sequenceid-USA/IA-CDC-LC1017569/2023|OQ631456.1|2023-02-16",
  //   name: "USA/IA-CDC-LC1017569/2023|OQ631456.1|2023-02-16",
  //   assemblyNames: ['NC_045512'],
  //   category: [],
  //   adapter: {
  //       type: "DashboardAdapter",
  //       bamLocation: {
  //           uri: "http://localhost:8000/uploads/my_vcf_haplotypes.bam"
  //       },
  //       index: {
  //           location: {
  //               uri: "http://localhost:8000/uploads/my_vcf_haplotypes.bam.bai"
  //           },
  //       },
  //   },
  //   displays: [
  //       {
  //           type:'LinearDashboardDisplay',
  //           displayId: 'sequence-display-id'+"USA/IA-CDC-LC1017569/2023|OQ631456.1|2023-02-16",
  //           groupname_tag: { 
  //             "group77": [{"unseenKey":"unaccounted1" , "mutation": "28877T:0.140006"}],
  //           }
  //       }
  //   ]
  // }

  //   {
  //     type: 'AlignmentsTrack',
  //     trackId: 'genes',
  //     name: 'spike-in_bams_file_0.bam',
  //     assemblyNames: ['NC_045512'],
  //     category: ['Genes'],
  //     adapter: {
  //       type: 'BamAdapter',
  //       bamLocation: {
  //         uri: 'spike-in_bams_file_0.bam',
  //       },
  //       index: {
  //         location: {
  //           uri: 'spike-in_bams_file_0.bam.bai',
  //         },
  //       },
  //     }
  // }
  ]
  
  export default tracks