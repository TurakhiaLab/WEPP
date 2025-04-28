import { lazy } from 'react'
import PluginManager from '@jbrowse/core/PluginManager'
import Plugin from '@jbrowse/core/Plugin'
import { ConfigurationSchema } from '@jbrowse/core/configuration'
import { DisplayType, TrackType, createBaseTrackConfig, createBaseTrackModel, RendererType } from '@jbrowse/core/pluggableElementTypes'

import LinearDashboardDisplayF from './LinearDashboardDisplay'

import AdapterType from '@jbrowse/core/pluggableElementTypes/AdapterType'
// import  DashboardRenderingSchema from './DashboardRenderer/configSchema'
// import DashboardRenderingComponent from './DashboardRenderer/component/DashboardRendering'
import WidgetType from '@jbrowse/core/pluggableElementTypes/WidgetType'
import { AdapterClass, configSchema as adapterSchema } from './DashboardAdapter'

import {DashboardRender, configSchema as DashboardRenderingSchema, ReactComponent as DashboardRenderingComponent} from './DashboardRenderer'

import {stateModelFactory as DashboardWidgetModelFactory, configSchema as dashboardWidgetConfigSchema } from './DashboardFeatureWidget'


export default class DashboardPlugin extends Plugin {
  name = 'DashboardPlugin'

  install(pluginManager: PluginManager) {
    const LGVPlugin = pluginManager.getPlugin(
      'LinearGenomeViewPlugin', 
    ) as import('@jbrowse/plugin-linear-genome-view').default

  const { BaseLinearDisplayComponent } = LGVPlugin.exports

  pluginManager.addTrackType(() => {
    const configSchema = ConfigurationSchema(
      'DashboardTrack',
      {}, 
      {
        baseConfiguration: createBaseTrackConfig(pluginManager),
        explicitIdentifier: 'trackId',
      }
    )
    const track = new TrackType({
      name: 'DashboardTrack',
      configSchema,
      stateModel: createBaseTrackModel(
        pluginManager, 
        'DashboardTrack',
        configSchema,
      ),
    })
    const linearPileupDisplay = pluginManager.getDisplayType('LinearPileupDisplay')
    const linearAlignmentDisplay = pluginManager.getDisplayType('LinearAlignmentsDisplay')
    const linearDashboardDisplay = pluginManager.getDisplayType('LinearDashboardDisplay')

    track.addDisplayType(linearPileupDisplay)
    track.addDisplayType(linearDashboardDisplay)
    track.addDisplayType(linearAlignmentDisplay)
    return track
  })

  pluginManager.addDisplayType(() => {
    const { configSchema, stateModel } = LinearDashboardDisplayF(pluginManager)
    return new DisplayType({
      name: 'LinearDashboardDisplay',
      configSchema,
      stateModel,
      trackType: 'DashboardTrack',
      viewType: 'LinearGenomeView',
      ReactComponent: BaseLinearDisplayComponent
    })
  })

  pluginManager.addAdapterType(
    () =>
      new AdapterType({
        name: 'DashboardAdapter',
        AdapterClass,
        configSchema: adapterSchema,
      }),
  )

  pluginManager.addRendererType(() => {
    return new DashboardRender({
      name: 'DashboardRenderer', 
      displayName: 'Dashboard renderer',
      ReactComponent: DashboardRenderingComponent, 
      configSchema: DashboardRenderingSchema, 
      pluginManager
    })
  })

    pluginManager.addWidgetType(
      () =>
        new WidgetType({
          name: 'DashboardFeatureWidget',
          heading: 'Feature details',
          configSchema: dashboardWidgetConfigSchema,
          stateModel: DashboardWidgetModelFactory(pluginManager),
          ReactComponent: lazy(() => import('./DashboardFeatureWidget/DashboardFeatureDetail')),
        }),
    )

  }
}