import { useCallback, useMemo, useEffect, useState } from "react";

// test
console.log("new worker");
//const workerPath = "../webworkers/localBackendWorker.js";

import workerSpec from "../webworkers/localBackendWorker.js?worker&inline";

//const url = new URL('../webworkers/localBackendWorker.js', import.meta.url)
//const getWorker = () => new Worker(url, { type: 'module' })

const worker = new workerSpec();

let onQueryReceipt = (receivedData) => {};
let onStatusReceipt = (receivedData) => {
  console.log("STATUS:", receivedData.data);
};

let onConfigReceipt = (receivedData) => {};
let onDetailsReceipt = (receivedData) => {};
let onListReceipt = (receivedData) => {};
let onNextStrainReceipt = (receivedData) => {
  console.log("NEXT STRAIN:", receivedData);
  // create a blob with this data and trigger download
  const blob = new Blob([JSON.stringify(receivedData)], {
    type: "application/json",
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.style.display = "none";
  a.href = url;
  a.download = "nextstrain.json";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  window.URL.revokeObjectURL(url);
};

let searchSetters = {};

worker.onmessage = (event) => {
  console.log(
    "got message from worker" //, event.data
  );
  if (event.data.type === "status") {
    onStatusReceipt(event.data);
  }
  if (event.data.type === "query") {
    onQueryReceipt(event.data.data);
  }
  if (event.data.type === "search") {
    // console.log("SEARCHRES", event.data.data);
    searchSetters[event.data.data.key](event.data.data);
  }
  if (event.data.type === "config") {
    onConfigReceipt(event.data.data);
  }
  if (event.data.type === "details") {
    onDetailsReceipt(event.data.data);
  }
  if (event.data.type === "list") {
    onListReceipt(event.data.data);
  }
  if (event.data.type === "nextstrain") {
    onNextStrainReceipt(event.data.data);
  }
};

function useLocalBackend(uploaded_data) {
  const [statusMessage, setStatusMessage] = useState({ message: null });
  onStatusReceipt = (receivedData) => {
    // console.log("STATUS:", receivedData.data);
    if (receivedData.data.error) {
      window.alert(receivedData.data.error);
      // console.log("ERROR33:", receivedData.data.error);
    }
    const total_nodes = receivedData.data.total;
    if (total_nodes && total_nodes > 6000000) {
      if (1) {
        window.alert(
          "This is a large tree which may use too much memory to run in the web browser. If the page crashes you might want to try the Taxonium desktop app."
        );
      }
    }
    setStatusMessage(receivedData.data);
  };
  useEffect(() => {
    console.log("Sending data to worker");
    worker.postMessage({
      type: "upload",
      data: uploaded_data,
    });
  }, [uploaded_data]);

  const queryNodes = useCallback(
    async (boundsForQueries, setResult, setTriggerRefresh, config) => {

      worker.postMessage({
        type: "query",
        bounds: boundsForQueries,
      });
      onQueryReceipt = (receivedData) => {
        //  console.log("CONFIG IS", config);

        receivedData.nodes.forEach((node) => {
          if (node.node_id === config.rootId) {
            node.mutations = config.rootMutations.map(
              (x) => config.mutations[x]
            );
          } else {
            node.mutations = node.mutations.map(
              (mutation) => config.mutations[mutation]
            );
          }
        });
        setResult(receivedData);
      };
    },
    []
  );

  const singleSearch = useCallback(
    (singleSearch, boundsForQueries, setResult) => {
      const key = JSON.parse(singleSearch).key;
      worker.postMessage({
        type: "search",
        search: singleSearch,
        bounds: boundsForQueries,
      });

      searchSetters[key] = (receivedData) => {

        setResult(receivedData);
      };
      return {
        abortController: {
          abort: () => console.log("no controller for local"),
        },
      };
    },
    []
  );

  const getDetails = useCallback((node_id, setResult) => {

    worker.postMessage({
      type: "details",
      node_id: node_id,
    });
    onDetailsReceipt = (receivedData) => {

      setResult(receivedData);
    };
  }, []);

  const getConfig = useCallback((setResult) => {

    worker.postMessage({
      type: "config",
    });

    onConfigReceipt = (receivedData) => {

      setResult(receivedData);
    };
  }, []);

  const getTipAtts = useCallback((nodeId, selectedKey, callback) => {
    
    worker.postMessage({
      type: "list",
      node_id: nodeId,
      key: selectedKey,
    });

    onListReceipt = (receivedData) => {

      callback(null, receivedData);
    };
  }, []);

  const getNextstrainJson = useCallback((nodeId, config) => {

    worker.postMessage({
      type: "nextstrain",
      node_id: nodeId,
      config: config,
    });
  }, []);

  return useMemo(() => {
    return {
      queryNodes,
      singleSearch,
      getDetails,
      getConfig,
      statusMessage,
      setStatusMessage,
      getTipAtts,
      getNextstrainJson,
      type: "local",
    };
  }, [
    queryNodes,
    singleSearch,
    getDetails,
    getConfig,
    statusMessage,
    setStatusMessage,
    getTipAtts,
    getNextstrainJson,
  ]);
}

export default useLocalBackend;
