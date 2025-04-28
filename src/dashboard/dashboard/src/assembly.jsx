import config from './config';


const assembly = {
    name: 'NC_045512',
    aliases: ['hg38'],
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'GRCh38-ReferenceSequenceTrack',
      adapter: {
        type: 'IndexedFastaAdapter',
        fastaLocation: {
          // uri: `${config.TAXONIUM_BASE}uploads/NC_045512v2.fa`,
          uri: `${config.REF_FA}`
        },
        faiLocation: {
          // uri: `${config.TAXONIUM_BASE}uploads/NC_045512v2.fa.fai`,
          uri: `${config.REF_FAI}`

        },
      },
    }
  }

  export default assembly
  