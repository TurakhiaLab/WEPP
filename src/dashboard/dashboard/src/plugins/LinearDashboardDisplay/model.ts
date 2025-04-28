import { cast, types } from 'mobx-state-tree'
import {
  AnyConfigurationSchemaType,
  ConfigurationReference,
  getConf
} from '@jbrowse/core/configuration'
import { getParentRenderProps } from '@jbrowse/core/util/tracks'
import {
  getSession,
  getContainingTrack,
  isSessionModelWithWidgets,
  Feature,
  getContainingView,
  getEnv,
  SimpleFeature,
} from '@jbrowse/core/util'
// icons
import FilterListIcon from '@mui/icons-material/FilterList'
import SortIcon from '@mui/icons-material/Sort'
import { LinearGenomeViewModel } from '@jbrowse/plugin-linear-genome-view'
import { SimpleFeatureSerialized } from '@jbrowse/core/util/simpleFeature'
import { getRpcSessionId } from '@jbrowse/core/util/tracks'

import { getUniqueModifications } from '../shared/getUniqueModifications'


import {
  createAutorun,
  getColorForModification,
  modificationData,
} from '../util'

import type {
  ModificationType,
  ModificationTypeWithColor,
  SortedBy,
} from '../shared/types'

// locals

import configSchemaF from './configSchema'
import PluginManager from '@jbrowse/core/PluginManager'
import {IFilter} from './FilterByTag'

import LinearGenomeViewPlugin from '@jbrowse/plugin-linear-genome-view'
import { observable } from 'mobx'
import { IAnyStateTreeNode, addDisposer, isAlive, getSnapshot } from 'mobx-state-tree'
import { lazy } from 'react'
import VisibilityIcon from '@mui/icons-material/Visibility'
import FilterAltIcon from '@mui/icons-material/FilterAlt';


const FilterByTagDialog = lazy(() => import('./FilterByTag'))
import { IAutorunOptions, autorun } from 'mobx'

type DisplayModel = IAnyStateTreeNode & { setError: (arg: unknown) => void }
import { getUniqueTagValues} from '../shared'
import { setFlagsFromString } from 'v8'
import { group } from 'console'
import { setDefaultHighWaterMark } from 'stream'
import { constrainedMemory } from 'process'
import { Checkbox } from '@headlessui/react'

  
  
  interface ExtraColorBy {
    custom?: Record<string, string>
  }
  
  const ColorByModel = types.maybe(
    types.model({
      type: types.string,
      tag: types.maybe(types.string),
      extra: types.frozen(),
    }),
  )
  const FilterModel = types.model({
    flagInclude: types.optional(types.number, 0),
    flagExclude: types.optional(types.number, 1540),
    readName: types.maybe(types.string),
    tagFilter: types.maybe(
      types.model({
        tag: types.string,
        value: types.string,
      }),
    ),
    filterReads: types.optional(types.frozen(), {})
  })

  const unseenModel = types.model({
    show: types.boolean,
    mutation: types.string
  })
  

  // using a map because it preserves order
const rendererTypes = new Map([
    ['pileup', 'PileupRenderer'],
    ['svg', 'SvgFeatureRenderer'],
    // ['Dashboard', 'DashboardRenderer']
  ])

type LGV = LinearGenomeViewModel

export default (
    pluginManager: PluginManager,
    configSchema: AnyConfigurationSchemaType,
  ) => {
    const { BaseLinearDisplay } = (
      pluginManager.getPlugin('LinearGenomeViewPlugin') as LinearGenomeViewPlugin
    )?.exports
  
    return types
      .compose(
        'LinearDashboardDisplay',
        BaseLinearDisplay,
        types.model({
          type: types.literal('LinearDashboardDisplay'),
          configuration: ConfigurationReference(configSchema),
          colorBy: ColorByModel,
          filterBy: types.optional(FilterModel, {}),
          showEPColor: false,
          unseenKeys: types.optional(types.map(unseenModel), {}),
          modificationsReady: false,
        }),
        
      )
      .volatile(() => ({
        colorTagMap: observable.map<string, string>({}),
        tagsReady: false,
        featureUnderMouseVolatile: undefined as undefined | Feature,
        group_id:'' as string,
        visibleModifications: observable.map<string, ModificationTypeWithColor>(
          {},
        ),
        
      }))
      .actions(self => ({
          setColorScheme(colorScheme: {
            type: string
            tag?: string
            // extra?: ExtraColorBy
          }) {
            self.colorTagMap = observable.map({}) // clear existing mapping
            self.colorBy = cast(colorScheme)
            if (colorScheme.tag) {
              self.tagsReady = false
            }
          },

          updateVisibleModifications(uniqueModifications: ModificationType[]) {
            uniqueModifications.forEach(value => {
              if (!self.visibleModifications.has(value.type)) {
                self.visibleModifications.set(value.type, {
                  ...value,
                  color: getColorForModification(value.type),
                })
              }
            })
          },

          setModificationsReady(flag: boolean) {
            self.modificationsReady = flag
          },

          updateColorTagMap(uniqueTag: string[]) {
            // pale color scheme
            // https://cran.r-project.org/web/packages/khroma/vignettes/tol.html
            // e.g. "tol_light"
            const colorPalette = [
              '#BBCCEE',
              'pink',
              '#CCDDAA',
              '#EEEEBB',
              '#FFCCCC',
              'lightblue',
              'lightgreen',
              'tan',
              '#CCEEFF',
              'lightsalmon',
            ]
    
            uniqueTag.forEach(value => {
              if (!self.colorTagMap.has(value)) {
                const totalKeys = [...self.colorTagMap.keys()].length
                self.colorTagMap.set(value, colorPalette[totalKeys])
              }
            })
          },

          setFilterBy(filter: IFilter) {
            self.filterBy = cast(filter)
          },

          setTagsReady(flag: boolean) {
            self.tagsReady = flag
          },
  
          setFeatureUnderMouse(feat?: Feature) {
            self.featureUnderMouseVolatile = feat
          },

          selectFeature(feature: Feature, display_id, all_group_name) {
            const session = getSession(self)
            if (isSessionModelWithWidgets(session)) {
              const featureWidget = session.addWidget(
                'DashboardFeatureWidget',
                'alignmentFeature',
                { featureData: feature.toJSON(), view: getContainingView(self),
                  track: getContainingTrack(self), 
                  "display_id": display_id,
                  'all_group_name': all_group_name
                 },
              )
              session.showWidget(featureWidget)
            }
            session.setSelection(feature)
          }
        }))

        .views(self => ({
          get autorunReady() {
            const view = getContainingView(self) as LGV
            return (
              view.initialized &&
              self.featureDensityStatsReady &&
              !self.regionTooLarge
            )
          },
        }))

        .actions(self => ({ 
          toggleEPDisplay() {
            self.showEPColor = !self.showEPColor
            if(self.showEPColor)
            {
              self.setColorScheme({type: 'tag', tag: 'EPP'})
            }
            else {
              self.setColorScheme({type: 'mappingQuality'});
            }
          },

          toggleUnaccountedMutationsDisplay(um_group: string) {

            
            const item = self.unseenKeys.get(um_group);
            if(item){
              item.show = !item.show;
            }

            const filter_reads_by_tag = {}

            filter_reads_by_tag[self.group_id] = []

            self.unseenKeys.forEach((value, key) => {
              if(value.show) {
                filter_reads_by_tag[self.group_id].push({'unseenKey':key, 'mutation': value.mutation})
              }
            });

            let um_filter_reads: IFilter = {
              flagExclude: 1540, 
              flagInclude: 0,
              filterReads: filter_reads_by_tag
            }
            self.setFilterBy(um_filter_reads)
        },

        afterAttach() {

          const is_sequence = self.configuration.toJSON().displayId.includes("sequence") ? true : false
          
          const group_name = getConf(self, 'groupname_tag')
          

            const group_name_keys = Object.keys(group_name);
            self.group_id = group_name_keys[0]



            group_name[self.group_id].map(item => {
              self.unseenKeys.set(item.unseenKey, unseenModel.create({ show: false, mutation: item.mutation}))
            })

            const filter_reads_tags = {}
            if(is_sequence){

              filter_reads_tags[self.group_id] = group_name[self.group_id]
            }
            else {
              filter_reads_tags[self.group_id] = []
            }
            

            self.setColorScheme({type: 'mappingQuality'});
            
            let filter_reads: IFilter = {
              flagExclude: 1540,
              flagInclude: 0,
              filterReads: filter_reads_tags
            }
            
            self.setFilterBy(filter_reads)

          createAutorun(
            self,
            async () => {

              const view = getContainingView(self) as LGV
              if (!self.autorunReady) {
                return
              }
              const { colorBy, tagsReady } = self

              const { staticBlocks } = view
              const vals = await getUniqueModifications({
                self,
                adapterConfig: getConf(self.parentTrack, 'adapter'),
                blocks: staticBlocks,
              })
              if (isAlive(self)) {
                self.updateVisibleModifications(vals)
                self.setModificationsReady(true)
              }

              if (colorBy?.tag && !tagsReady) {
                const vals = await getUniqueTagValues(self, colorBy, staticBlocks)
                self.updateColorTagMap(vals)
                self.setTagsReady(true)
              }
              self.setTagsReady(true)
            },
            { delay: 1000 },
          )

          addDisposer(
            self,
            autorun(async () => {
              const session = getSession(self)
              try {
                const featureId = self.featureIdUnderMouse
                if (self.featureUnderMouse?.id() !== featureId) {
                  if (!featureId) {
                    self.setFeatureUnderMouse(undefined)
                  } else {
                    const sessionId = getRpcSessionId(self)
                    const view = getContainingView(self)
                    const { feature } = (await session.rpcManager.call(
                      sessionId,
                      'CoreGetFeatureDetails',
                      {
                        featureId,
                        sessionId,
                        layoutId: view.id,
                        rendererType: 'PileupRenderer',
                      },
                    )) as { feature: SimpleFeatureSerialized }
  
                    // check featureIdUnderMouse is still the same as the
                    // feature.id that was returned e.g. that the user hasn't
                    // moused over to a new position during the async operation
                    // above
                    if (self.featureIdUnderMouse === feature?.uniqueId) {
                      self.setFeatureUnderMouse(new SimpleFeature(feature))
                    }
                  }
                }
              } catch (e) {
                console.error(e)
                session.notify(`${e}`);
              }
            }),
          )
      }

    }
    ))
      .views(self => ({

        get visibleModificationTypes() {
          return [...self.visibleModifications.keys()]
        },

        get rendererConfig() {
          const { rendererTypeName } = self
          const configBlob = getConf(self, ['renderers', rendererTypeName]) || {}
          return self.rendererType.configSchema.create(
            {
              ...configBlob,
            },
            getEnv(self),
          )
        },
        renderReady() {
          return self.tagsReady
        },
        })
        )
      .views(self => {
          const {
            renderProps: superRenderProps,
            trackMenuItems: superTrackMenuItems,
          } = self
          return {
            renderProps() {
              return {
                ...superRenderProps(),
                ...getParentRenderProps(self),
                displayModel: self,
                config: self.rendererConfig,
              }
            },
  
          get rendererTypeName() {
              const viewName = getConf(self, 'defaultRendering')
              const rendererType = rendererTypes.get(viewName)
              if (!rendererType) {
                throw new Error(`unknown alignments view name ${viewName}`)
              }
              return rendererType
            },

          trackMenuItems() {
            return [
              ...superTrackMenuItems(),
              //   {
              //     label: 'Filter by',
              //     icon: FilterListIcon,
              //     onClick: () => {
              //       getSession(self).queueDialog(doneCallback => [
              //         FilterByTagDialog,
              //         { model: self, handleClose: doneCallback },
              //       ])
              //     },
              // },
              {
                label: 'color by EPP',
                icon: VisibilityIcon,
                type: 'checkbox',
                checked: self.showEPColor,
                onClick: () => {
                  self.toggleEPDisplay()
                }
              },
              {
                label: "unaccounted mutations",
                icon: SortIcon,
                subMenu: Array.from(self.unseenKeys.entries()).map(([key, item]) => ({
                  label: item.mutation.split(":")[0],
                  type:'checkbox',
                  checked: item.show,
                  onClick: () => {
                    self.toggleUnaccountedMutationsDisplay(key)
                  }
                }))
              }
            ]
            },

            renderPropsPre() {
              const { colorTagMap, colorBy, filterBy, rpcDriverName, visibleModifications } = self

              const all_groups = getConf(self, 'all_group_name')

              const superProps = superRenderProps()
              return {
                ...superProps,
                notReady: superProps.notReady || !self.renderReady(),
                rpcDriverName,
                displayModel: self,
                visibleModifications: Object.fromEntries(
                  visibleModifications.toJSON(),
                ),
                colorBy: colorBy ? getSnapshot(colorBy) : undefined,
                filterBy: JSON.parse(JSON.stringify(filterBy)),
                colorTagMap: Object.fromEntries(colorTagMap.toJSON()),
                config: self.rendererConfig,
                async onFeatureClick(_: unknown, featureId?: string) {
                  const session = getSession(self)
                  const { rpcManager } = session
                  try {
                    const f = featureId || self.featureIdUnderMouse
                    if (!f) {
                      self.clearFeatureSelection()
                    } else {
                      const sessionId = getRpcSessionId(self)
                      const { feature } = (await rpcManager.call(
                        sessionId,
                        'CoreGetFeatureDetails',
                        {
                          featureId: f,
                          sessionId,
                          layoutId: getContainingView(self).id,
                          rendererType: 'PileupRenderer',
                        },
                      )) as { feature: SimpleFeatureSerialized | undefined }
    
                      if (feature) {

                        self.selectFeature(new SimpleFeature(feature), self.configuration.toJSON().displayId, all_groups)
                      }
                    }
                  } catch (e) {
                    console.error(e)
                    session.notify(`${e}`)
                  }
                },
                onClick() {
                  self.clearFeatureSelection()
                },
                // similar to click but opens a menu with further options
                async onFeatureContextMenu(_: unknown, featureId?: string) {
                  const session = getSession(self)
                  const { rpcManager } = session
                  try {
                    const f = featureId || self.featureIdUnderMouse
                    if (!f) {
                      self.clearFeatureSelection()
                    } else {
                      const sessionId = getRpcSessionId(self)
                      const { feature } = (await rpcManager.call(
                        sessionId,
                        'CoreGetFeatureDetails',
                        {
                          featureId: f,
                          sessionId,
                          layoutId: getContainingView(self).id,
                          rendererType: 'PileupRenderer',
                        },
                      )) as { feature: SimpleFeatureSerialized }
    
                      if (feature) {
                        self.setContextMenuFeature(new SimpleFeature(feature))
                      }
                    }
                  } catch (e) {
                    console.error(e)
                    session.notify(`${e}`)
                  }
                }
              }
            }
          }
        })
        .views(self => ({
          renderProps() {
            return self.renderPropsPre()
          },
        }))
  }
  
