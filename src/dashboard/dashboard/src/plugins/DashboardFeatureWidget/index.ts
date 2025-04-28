import { lazy } from 'react'

import WidgetType from '@jbrowse/core/pluggableElementTypes/WidgetType'

import { ConfigurationSchema } from '@jbrowse/core/configuration'

export const configSchema = ConfigurationSchema('DashboardFeatureWidget', {})

// import { configSchema } from './configSchema'
export { stateModelFactory } from './stateModelFactory'
