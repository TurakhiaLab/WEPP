import React, { useState, useEffect, useMemo, useCallback, useRef } from 'react'

import { getDefaultSearch } from "./utils/searchUtil.js";

import Taxonium from './Taxonium.jsx'
import tracks from './tracks.jsx';
import assembly from './assembly.jsx';
import defaultSession from './defaultSession.jsx'
import { sourceData } from './data.js';
// import "@jbrowse/plugin-alignments"   

import {
    createViewState,
    JBrowseLinearGenomeView,
  } from '@jbrowse/react-linear-genome-view'

  const createDefaultSearch = (mark_nodes) => {
    return mark_nodes.map(node => {
      return getDefaultSearch(null, node);
    })
  }

  const App = () => {

    const [viewState, setViewState] = useState();
    const [all_tracks, setTracks] = useState([])
    const [patches, setPatches] = useState();
    const [showTrack, setShowTrack] = useState(null)

    const clickedNodeRef = useRef(null);
    const queryRef = useRef(null);
    const mark_nodeRef = useRef(null);
    const trackIDsRef = useRef([]);
    
    const [createTract, setCreateTrack] = useState(null)

    const [JBrowseOpen, setJBrowseOpen] = useState(false);
    const [selectedFile, setSelectedFile] = useState(true);

    useEffect(() => {

        const state = createViewState({
          assembly,
          tracks: tracks,
          defaultSession,
          onChange: patch => {
            setPatches(previous => previous + JSON.stringify(patch) + '\n')
          },
        })
        setViewState(state)
      }, []);
      
    const onClickNode = useCallback((selectedNode) => {
    if ( selectedNode && mark_nodeRef.current.includes(selectedNode.nodeDetails.name)){
            // addTrack(showTrack.trackID, showTrack.bam_location, showTrack.bami_location)
            clickedNodeRef.current = selectedNode
            
            selectedNode.clearNodeDetails()
            // addTrack();
            
            setCreateTrack(true);

            setJBrowseOpen(true);
    }
    }, [JBrowseOpen]);


    const addTrack = () => {

    const trackId = clickedNodeRef.current.nodeDetails.name +"sdaasd"
    const bam_location = "/data/"+selectedFile[clickedNodeRef.current.nodeDetails.name]
    const bami_location = bam_location+".bai";

    const new_track_addition = {
        type: "AlignmentsTrack",
        trackId: trackId,
        name: 'bamfile-name-'+clickedNodeRef.current.nodeDetails.name,
        assemblyNames: ['NC_045512'],
        category: [],
        adapter: {
        type: "BamAdapter",
        bamLocation: {
            uri: bam_location
        },
        index: {
            location: {
            uri: bami_location
            },
        },
        }
    }
    
    if (!trackIDsRef.current.includes(trackId))
    {
        setTracks([...all_tracks, new_track_addition])
        viewState.session.addTrackConf(new_track_addition);
    
        trackIDsRef.current = [...trackIDsRef.current, trackId];
        setShowTrack(trackId);
    }

    }

useEffect(() => {
    if(createTract){
        addTrack();
    }      
    }, [createTract])

    useEffect(() => {
    if(showTrack)
    {
        trackIDsRef.current.forEach((trackId) => {
        viewState.session.view.showTrack(trackId);
        })
        
        setShowTrack(null)
        setCreateTrack(false);
    }
    }, [showTrack])

    const handleFileUpload = (event) => {
        const array = ["ASV_2035", "ASV_1712", "ASV_2141", "ASV_1647"];
  
        const fileNames = Array.from(event.target.files).map(file => file.name);
        
        const mark_nodes_s = array.sort(() => 0.5 - Math.random()).slice(0, fileNames.length);
        mark_nodeRef.current = mark_nodes_s;

        const nodeFileDictionary = {};
        mark_nodes_s.forEach((node, index) => {
          if (fileNames[index]) {
            nodeFileDictionary[node] = fileNames[index];
          }
        });
        
        const default_search = createDefaultSearch(mark_nodes_s);
  
        const zoom_to_indexes = []
        for (let i = 0; i < fileNames.length; i++) {
          zoom_to_indexes.push(i.toString());
        }

        queryRef.current = {
          srch: JSON.stringify(default_search),
          enabled: JSON.stringify(Object.fromEntries(default_search.map(value => [value.key, true]))),
          backend: "",
          xType: "x_dist",
          zoomToSearch: zoom_to_indexes,
          mutationTypesEnabled: JSON.stringify({ aa: true, nt: false }),
          treenomeEnabled: false,
        };        
        setSelectedFile(nodeFileDictionary);
      }
  
if (!viewState) {
    return null
}
  
const handleUpload = (event) => {
    console.log("File Uploaded");
  }


return (
    <>
    <div>
          <div style={{ display: "flex", justifyContent: "center" }}>
              <h1 style={{ textAlign: "center", fontSize: "24px", fontWeight: "bold", color: "#333" }}>Waste Water Dashboard</h1>
          </div>
          <br />
            {!selectedFile ? 
                <div className='flex' style={{ justifyContent: "center", alignItems: "center", height: "100vh" }}>
                <input type="file" onChange={handleFileUpload} style={{ padding: '10px', margin: '5px', border: '1px solid #ccc', borderRadius: '5px' }} multiple />
                <button onClick={handleUpload} style={{ padding: '10px 20px', margin: '5px', backgroundColor: '#4CAF50', color: 'white', border: 'none', borderRadius: '5px', cursor: 'pointer' }}>Upload</button>
                </div>
          :
         <div>
        <div style={{display:"flex", height: "80vh"}} >
          <div className='flex-grow'>
            <Taxonium sourceData={sourceData} query={queryRef.current} onClickNode={onClickNode}/>
            {/* <Taxonium backendUrl="https://api.cov2tree.org" query={queryRef.current} onClickNode={onClickNode} />  */}

          </div>
          <div className='flex' style={{margin: "10px"}}>
            {JBrowseOpen &&
              <div>
                {/* <button onClick={() => {
                  setJBrowseOpen(false) 
                  clickedNodeRef.current.clearNodeDetails()
                  }}>X</button> */}
                <JBrowseLinearGenomeView viewState={viewState} />
              </div>
          }
          </div> 
        </div> 
        </div> 
      }
        </div>
    </>
)
}

  export default App