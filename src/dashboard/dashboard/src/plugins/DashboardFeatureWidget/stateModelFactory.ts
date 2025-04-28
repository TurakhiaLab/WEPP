import { stateModelFactory as baseModelFactory } from '@jbrowse/core/BaseFeatureWidget'
import { types } from 'mobx-state-tree'

import type PluginManager from '@jbrowse/core/PluginManager'
import type { Instance } from 'mobx-state-tree'

export function stateModelFactory(pluginManager: PluginManager) {
  const baseModel = baseModelFactory(pluginManager)
  return types.compose(
    baseModel,
    
    types.model('DashboardFeatureWidget', {
      type: types.literal('DashboardFeatureWidget'),
      display_id: types.string,
      all_group_name: types.frozen(),
    }),
  )
}

export type AlignmentFeatureWidgetStateModel = ReturnType<
  typeof stateModelFactory
>

export type AlignmentFeatureWidgetModel =
  Instance<AlignmentFeatureWidgetStateModel>