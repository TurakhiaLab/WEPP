const assembly = {
    name: 'NC_045512',
    aliases: ['hg38'],
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'GRCh38-ReferenceSequenceTrack',
      adapter: {
        type: 'IndexedFastaAdapter',
        fastaLocation: {
          uri: '/konstantinos/refsequences.fasta',
        },
        faiLocation: {
          uri: '/konstantinos/refsequences.fasta.fai',
        },
      },
    }
  }

  export default assembly
  