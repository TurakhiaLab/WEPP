import PluginManager from '@jbrowse/core/PluginManager'
import { ConfigurationSchema } from '@jbrowse/core/configuration'
import LinearGenomeViewPlugin from '@jbrowse/plugin-linear-genome-view'
import {types} from 'mobx-state-tree'

export default (pluginManager: PluginManager) => {

    const { baseLinearDisplayConfigSchema } = (
    pluginManager.getPlugin('LinearGenomeViewPlugin') as LinearGenomeViewPlugin
  ).exports
  
  // Define the unseen object schema
  const Unseen = types.model({
    unseenKey: types.string,    // e.g., "unseen1" or "unseen2"
    value: types.string         // e.g., "AT10:10%" or "GC10:20%"
  });

  // Define the group schema which holds an array of unseen objects
  const Group = types.array(Unseen);

  // Define the root schema
  const RootModel = types.model({
    groupname_tag: types.map(types.array(Unseen))    // Map where the keys are group names and values are arrays of unseen objects
  });

  return ConfigurationSchema(
    'LinearDashboardDisplay',
    {
      defaultRendering: {
        type: 'stringEnum',
        model: types.enumeration('Rendering', ['pileup', 'Dashboard']),
        defaultValue: 'pileup',
    },
    renderers: ConfigurationSchema('RenderersConfiguration', {
      DashboardRenderer:
      pluginManager.getRendererType('DashboardRenderer').configSchema,
    }),

    colorScheme: {
      type: 'stringEnum',
      model: types.enumeration('colorScheme', [
        'strand',
        'normal',
        'insertSize',
        'insertSizeAndOrientation',
        'mappingQuality',
        'tag',
      ]),
      description: 'color scheme to use',
      defaultValue: 'mappingQuality',
    },
    all_group_name: {
      type: 'frozen',
      defaultValue: {}
    },
    groupname_tag: {
      // type: 'string',
      // defaultValue: "",
      type: 'frozen',
      mode: RootModel,
      description: "read groups",
      defaultValue: {
        // "group1": [
        //   {unseenKey:'unseen8',mutation:'22599A:0.18'},
        //   {unseenKey:'unseen5',mutation:'21987A:0.55'},
        //   {unseenKey:'unseen4',mutation:'26364C:0.10'},
        //   {unseenKey:'unseen3',mutation:'1452A:0.20'},
        // ],
      },
    },
  },
    { baseConfiguration: baseLinearDisplayConfigSchema, explicitlyTyped: true },
  )
}