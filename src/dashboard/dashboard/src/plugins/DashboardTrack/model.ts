
import {
    AnyConfigurationSchemaType
  } from '@jbrowse/core/configuration'

import { types } from 'mobx-state-tree'

import { createBaseTrackModel } from '@jbrowse/core/pluggableElementTypes'
import PluginManager from '@jbrowse/core/PluginManager'

export default (
    pluginManager: PluginManager,
    baseTrackConfig: AnyConfigurationSchemaType,
  ) => {
    return types.
    compose(
        'DashboardTrack',
        createBaseTrackModel(
            pluginManager, 
            'DashboardTrack',
            baseTrackConfig,
          ),
          types.model({
            groupname: types.optional(types.string, "", [null]),
          }),
    )
  }