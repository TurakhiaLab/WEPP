import { ConfigurationSchema } from '@jbrowse/core/configuration'
import { types } from 'mobx-state-tree'


export default ConfigurationSchema(
  'DashboardAdapter',
  {
    bamLocation: {
        type: 'fileLocation',
        defaultValue: { uri: '/path/to/my.bam', locationType: 'UriLocation' },
      },
    index: ConfigurationSchema('BamIndex', {
        indexType: {
            model: types.enumeration('IndexType', ['BAI', 'CSI']),
            type: 'stringEnum',
            defaultValue: 'BAI',
            },

        location: {
            type: 'fileLocation',
            defaultValue: {
                uri: '/path/to/my.bam.bai',
                locationType: 'UriLocation',
            },
        },
    })
},
{ explicitlyTyped: true },
)