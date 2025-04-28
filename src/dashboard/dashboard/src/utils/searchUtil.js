export function getDefaultSearch(config, text, hs_value ) {
    const key = Math.random().toString(36).substring(2, 15);
    // if (!key) {
    //   key = Math.random().toString(36).substring(2, 15);
    //   console.log("generated key", key);
    // }
    if (config && config.defaultSearch) {
      return config.defaultSearch;
    }
    return {
      key,
      type: "name",
      method: "text_exact",
      text,
      gene: "S",
      position: 484,
      new_residue: "any",
      min_tips: 0,
      ...(hs_value ? { hs_value: (parseFloat(hs_value) * 100).toFixed(3) } : {})
    };
  }
  