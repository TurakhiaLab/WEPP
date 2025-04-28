import React from 'react'
import { Paper } from '@mui/material'
import { observer } from 'mobx-react'
import clone from 'clone'
import { FeatureDetails } from '@jbrowse/core/BaseFeatureWidget/BaseFeatureDetail'

// locals
import { getTag } from './util'
import { tags } from './taginfo'
import { AlignmentFeatureWidgetModel } from './stateModelFactory'

// local components
import SuppAlignments from './SuppAlignments'
import UnnaccountedMutations from './UnaccountedMutationsTable'
import Haplotypes from './haplotypes'
import Flags from './Flags'
import PairLink from './PairLink'
import Formatter from './Formatter'
import {
    SimpleFeature,
  } from '@jbrowse/core/util'

const omit = ['clipPos', 'flags', "UM", "RG"]

const AlignmentsFeatureDetails = observer(function (props: {
  model: AlignmentFeatureWidgetModel
}) {
  const { model } = props
  const feat = clone(model.featureData)

  const display_id = model.display_id
  const all_group_name = model.all_group_name
  if (display_id.includes("sequence")) {
    omit.push('qual', 'template_length', 'Score', 'MQ', 'CIGAR', 'length_on_ref', 'seq_length', 'Type')
  }

  const SA = getTag('SA', feat) as string

  const UM = getTag('UM', feat) as string
  let RG = getTag('RG', feat) as string

  const RG_arrays = RG?.split(",")
  
  return (
    <Paper data-testid="alignment-side-drawer">
      <FeatureDetails
        {...props}
        omit={omit}
        // @ts-expect-error
        descriptions={{ ...tags, tags: tags }}
        feature={feat}
        formatter={(value, key) =>
          key === 'next_segment_position' ? (
            <PairLink model={model} locString={value as string} />
          ) : (
            <Formatter value={value} />
          )
        }
      />
      
      {!display_id.includes("sequence") ? (
        <Haplotypes haplotype_names={all_group_name} RG={RG_arrays} />
      ): (
        <>/</>
      )
    }

      <UnnaccountedMutations UM={UM} is_sequence={display_id.includes("sequence")}/>

      {SA ? <SuppAlignments model={model} tag={SA} feature={feat} /> : null}
      {feat.flags !== undefined ? <Flags feature={feat} {...props} /> : null}
    </Paper>
  )
})

export default AlignmentsFeatureDetails