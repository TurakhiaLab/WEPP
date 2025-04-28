import config from '../config';

const addTrack = (clickedNodeRef, selectedFile, trackIDsRef, setTracks, setShowTrack, viewState, all_tracks) => {

    // {'group1': {'filename': 'output_multi.bam', 'node_name': 'Germany/IMS-10209-CVDP-D48209F5-5BED-436E-BFC4-D2118C232BC4/2021'}, 
    // 'group2': {'filename': 'output_multi.bam', 'node_name': 'England/PHEC-Z306ZA27/2021'}}
    const node_name = clickedNodeRef.current.nodeDetails.name

    // const key = Object.keys(selectedFile).find(key => selectedFile[key].node_name === node_name);

    const read_groupname = selectedFile[node_name]['groupname']

    const unseenKey_dict = []

    const all_group_name = {}
    Object.entries(selectedFile).forEach(([node_name, read_groups]) => {
        all_group_name[read_groups['groupname']] = node_name
    })

    selectedFile[node_name][read_groupname].forEach(obj => {
        Object.entries(obj).forEach(([key, value]) => {
            unseenKey_dict.push({
                "unseenKey":key, 
                "mutation":value
            })
        });
    });
    
    // const read_groupname = key

    const bam_filename = selectedFile[node_name]['filename']
    const trackId = node_name

    // const bam_location = `${config.TAXONIUM_BASE}/uploads/${bam_filename}`
    const bam_location = `${config.BAM}${bam_filename}`
    const bami_location = bam_location + ".bai";
    
    const new_track_addition = {
        type: "DashboardTrack",
        trackId: trackId,
        name: "Reads mapping to "+clickedNodeRef.current.nodeDetails.name,
        assemblyNames: ['NC_045512'],
        category: [],
        adapter: {
            type: "DashboardAdapter",
            bamLocation: {
                uri: bam_location
            },
            index: {
                location: {
                    uri: bami_location
                },
            },
        },
        displays: [
            {
                type:'LinearDashboardDisplay',
                displayId: 'display-id'+trackId,
                groupname_tag: { 
                    [read_groupname]: unseenKey_dict
                },
                all_group_name : all_group_name
                // groupname_tag: read_groupname
                /*
                    groupname_tag: {
                        group1: [{unseen1: "AT10:10%"}],
                        group2: [{unseen2: "GC10:20%"}]
                        }
                    }
                */
            }
        ]
    };


    const haplotype_sequences = selectedFile['HP_SEQ']['filename']



    // const hap_location = `${config.TAXONIUM_BASE}/uploads/${haplotype_sequences}`;
    const hap_location = `${config.BAM}${haplotype_sequences}`;
    
    const hapi_location = hap_location + ".bai";

    const haplotype_trackid = "sequenceid-"+read_groupname

    const haplotype_track = {

        type: "DashboardTrack",
        trackId: haplotype_trackid,
        name: "sequence-"+clickedNodeRef.current.nodeDetails.name,
        assemblyNames: ['NC_045512'],
        category: [],
        adapter: {
            type: "DashboardAdapter",
            bamLocation: {
                uri: hap_location
            },
            index: {
                location: {
                    uri: hapi_location
                },
            },
        },
        displays: [
            {
                type:'LinearDashboardDisplay',
                displayId: 'sequence-display-id'+"haplotype_track"+read_groupname,
                groupname_tag: {
                    [read_groupname]: unseenKey_dict
                },
                all_group_name : all_group_name
            }
        ]
    }

    if (!trackIDsRef.current.includes(trackId))
    {
        // setTracks([...all_tracks, new_track_addition, haplotype_track])
        // setTracks([...all_tracks, new_track_addition])
        viewState.session.addTrackConf(haplotype_track);
        viewState.session.addTrackConf(new_track_addition);

        viewState.session.view.showTrack(haplotype_trackid);
        viewState.session.view.showTrack(trackId);

        trackIDsRef.current = [...trackIDsRef.current, haplotype_trackid, trackId];
        // trackIDsRef.current = [...trackIDsRef.current, trackId];
        // setShowTrack(trackId);
    }
    else {
        viewState.session.view.showTrack(haplotype_trackid);
        viewState.session.view.showTrack(trackId);
    }
}

export default addTrack