
const DEV = false

let API_BASE = null
let TAXONIUM_BASE = null

if(DEV) {
  API_BASE = 'http://localhost:8000/';
  TAXONIUM_BASE = 'http://localhost:8080/';
}else {
  API_BASE = window.location.origin + '/api';
  TAXONIUM_BASE = window.location.origin + '/taxonium/'
}

const REF_FA = `${TAXONIUM_BASE}uploads/NC_045512v2.fa`
const REF_FAI = `${TAXONIUM_BASE}uploads/NC_045512v2.fa.fai`

const BAM = `${TAXONIUM_BASE}uploads/`

const PROJECTS = `${TAXONIUM_BASE}api/projects`

const RESULT = `${TAXONIUM_BASE}api/result/`

const UPLOAD = `${TAXONIUM_BASE}api/upload`

export default {
  API_BASE,
  TAXONIUM_BASE,
  REF_FA,
  REF_FAI,
  BAM,
  PROJECTS,
  RESULT,
  UPLOAD
};