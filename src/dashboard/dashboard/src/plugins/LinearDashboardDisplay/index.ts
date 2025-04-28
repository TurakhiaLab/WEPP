import PluginManager from "@jbrowse/core/PluginManager";
import configSchemaF from "./configSchema";
import modelF from './model'

export default (pluginManager: PluginManager) => {
    const schema = configSchemaF(pluginManager)
    return {
      configSchema: schema,
      stateModel: modelF(pluginManager, schema),
    }
  }
