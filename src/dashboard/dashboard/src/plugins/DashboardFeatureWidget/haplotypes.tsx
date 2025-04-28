import React from 'react'
import { BaseCard } from '@jbrowse/core/BaseFeatureWidget/BaseFeatureDetail'
import { SimpleFeatureSerialized, getEnv, getSession } from '@jbrowse/core/util'
import { ViewType } from '@jbrowse/core/pluggableElementTypes'

// locals
import { AlignmentFeatureWidgetModel } from './stateModelFactory'
import SuppAlignmentsLocStrings from './SuppAlignmentsLocStrings'
import LaunchBreakpointSplitViewPanel from './LaunchBreakpointSplitViewPanel'


export default function Haplotypes(props) {

  const {haplotype_names, RG } = props
  console.log("haplotype_names", haplotype_names, RG)
  return (
    <BaseCard title="Other Haplotypes">
    <div style={{ maxHeight: '400px', overflowY: 'auto' }}>
    <ul style={{ listStyleType: 'none', padding: '0' }}>
      {RG && RG.length>0 ? (RG.map((haplotype) => (
        <li key={haplotype} style={{ backgroundColor: '#f9f9f9', borderBottom: '1px solid #ddd', padding: '8px' }}>
          <span style={{ fontWeight: 'bold' }}>{haplotype_names[haplotype]}</span>
        </li>
      ))) : (
        <span style={{ fontWeight: 'bold' }}>No Haplotypes</span>
      )
    }
    </ul>
    </div>
    </BaseCard>
  )
}