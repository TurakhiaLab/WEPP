import React, { useState, Suspense, useRef, useEffect, useCallback } from "react";
import TaxoniumBit from "./components/TaxoniumBit";
import { getDefaultSearch } from "./utils/searchUtil.js";

import tracks from './tracks.jsx';
import assembly from './assembly.jsx';
import defaultSession from './defaultSession.jsx'
import addTrack  from './utils/TrackUtils.jsx'
import { Header } from './utils/UIUtils.jsx'
import FileUpload from './utils/uploadUtils.jsx'
import DashboardPlugin from './plugins'
import config from './config';

import {
  createViewState,
  JBrowseLinearGenomeView,
} from '@jbrowse/react-linear-genome-view'

const default_query = {};

default_query.backend = null;


const createDefaultSearch = (mark_nodes) => {
  
  return Object.keys(mark_nodes).map(node => {
    return getDefaultSearch(null, node, mark_nodes[node]);
  })
}

function App() {
  // Import taxonium-component
  useEffect(() => {
    import("taxonium-component");
  }, []);

  const [viewState, setViewState] = useState();
  const [all_tracks, setTracks] = useState([])
  const [showTrack, setShowTrack] = useState(null)
  const clickedNodeRef = useRef(null);
  const mark_nodeRef = useRef(null);
  const trackIDsRef = useRef([]);
  const [createTract, setCreateTrack] = useState(null)
  const [JBrowseOpen, setJBrowseOpen] = useState(false);
  const [selectedFile, setSelectedFile] = useState(false);
  const [uploadProgress, setuploadProgress] = useState(0);
  // const [query, updateQuery] = useState(default_query);
  
  const [patches, setPatches] = useState();
  
  const [backupQuery, setBackupQuery] = useState(default_query);
  const backupUpdateQuery = useCallback((newQuery) => {
    setBackupQuery((oldQuery) => ({ ...oldQuery, ...newQuery }));
  }, []);
  
  const query = backupQuery; // query
  const updateQuery = backupUpdateQuery;

  const [refNames, setRefNames] = useState([]);

  const toggleJBrowse = () => {
    setJBrowseOpen(!JBrowseOpen);
  };

  useEffect(() => {
    const state = createViewState({
      assembly,
      tracks: tracks,
      defaultSession,
      plugins:[
        DashboardPlugin
        ],
      onChange: patch => {
        setPatches(previous => previous + JSON.stringify(patch) + '\n')
      },
    })
    setViewState(state)
  }, []);

useEffect(() => {
  if(createTract){
      addTrack(clickedNodeRef, 
        selectedFile, 
        trackIDsRef, 
        setTracks, 
        setShowTrack, 
        viewState, all_tracks);
      }
    setCreateTrack(false);
  }, [createTract])

  useEffect(() => {
    if(showTrack)
    {
        // trackIDsRef.current.forEach((trackId) => {
        // viewState.session.view.showTrack(trackId);
        // })
        viewState.session.view.showTrack(showTrack);

        setShowTrack(null)
        setCreateTrack(false);
    }
    }, [showTrack])

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
    if (!viewState) {
      return null
  }

  // console.log("Taxonium", config.TAXONIUM_BASE)
  return (
    <>
    {/* <JBrowseLinearGenomeView viewState={viewState} /> */}
      <div>
          < Header/>
          <br />
            {!selectedFile ? 
              <div>
                <FileUpload setSelectedFile={setSelectedFile} createDefaultSearch={createDefaultSearch} mark_nodeRef={mark_nodeRef} updateQuery={updateQuery}/>
              </div>
              
          :
         <div>
        <div style={{display:"flex", height: "100vh"}} >
          <div className="h-screen w-screen flex flex-col overflow-hidden">
            <div className="h-full">
              <TaxoniumBit
                backendUrl={`${config.TAXONIUM_BASE}`}
                // backendUrl="http://localhost:8080"
                // backendUrl="https://api.cov2tree.org"
                query={query} updateQuery={updateQuery} onClickNode={onClickNode}
              />
            </div>
        </div>
          <div className='flex' style={{margin: "10px"}}>
            {JBrowseOpen &&
            <div>
              <JBrowseLinearGenomeView viewState={viewState} />
            </div>
            }
          </div>
        </div> 
      </div>
      }
    </div>
  </>
);
}
export default App;