class m extends Map {
  constructor(e, t = { checkIndent: !0 }) {
    super();
    const { checkIndent: n } = t;
    this._checkIndent = n;
    let i;
    typeof e == "string" ? i = e.trimEnd().split(/\r?\n/) : e ? i = e : i = [], this._keyAndCommentOrder = [], i.forEach((s) => {
      this.add(s);
    });
  }
  /**
   * Add a single line to the stanza. If the exact line already exists, does
   * nothing.
   * @param {string} line A stanza line
   * @returns {RaStanza} The RaStanza object
   */
  add(e) {
    if (e === "")
      throw new Error("Invalid stanza, contained blank lines");
    if (e.trim().startsWith("#"))
      return this._keyAndCommentOrder.push(e.trim()), this;
    if (e.trimEnd().endsWith("\\")) {
      const h = e.trimEnd().slice(0, -1);
      return this._continuedLine ? this._continuedLine += h.trimStart() : this._continuedLine = h, this;
    }
    let t = e;
    if (this._continuedLine && (t = this._continuedLine + t.trimStart(), this._continuedLine = void 0), this.indent || this._checkIndent) {
      const h = t.match(/^([ \t]+)/);
      if (this.indent === void 0)
        h ? [, this.indent] = h : this.indent = "";
      else if (this.indent === "" && h !== null || this.indent && h && this.indent !== h[1])
        throw new Error("Inconsistent indentation of stanza");
    } else
      this.indent = "";
    const n = t.trim(), i = n.indexOf(" ");
    if (i === -1) {
      if (!this.nameKey)
        throw new Error("First line in a stanza must have both a key and a value");
      return this.has(n) ? this : (this._keyAndCommentOrder.push(n), super.set(n, ""));
    }
    const s = n.slice(0, i), r = n.slice(i + 1);
    if (this.has(s) && r !== this.get(s))
      throw new Error(`Got duplicate key with a different value in stanza: "${s}" key has both ${this.get(s)} and ${r}`);
    return this._keyAndCommentOrder.push(s), this.nameKey || (this.nameKey = s, this.name = n.slice(i + 1)), super.set(s, r);
  }
  /**
   * Use `add()` if possible instead of this method. If using this, be aware
   * that no checks are made for comments, indentation, duplicate keys, etc.
   * @param {string} key The key of the stanza line
   * @param {string} value The value of the stanza line
   * @returns {RaStanza} The RaStanza object
   */
  set(e, t) {
    if (typeof t != "string")
      throw new Error(`Value of ${e} must be a string, got ${typeof t}`);
    return super.set(e, t);
  }
  /**
   * Delete a line
   * @param {string} key The key of the line to delete
   * @returns {boolean} true if the deleted line existed, false if it did not
   */
  delete(e) {
    if (e === this.nameKey)
      throw new Error("Cannot delete the first line in a stanza (you can still overwrite it with set()).");
    return this._keyAndCommentOrder.includes(e) && (this._keyAndCommentOrder = this._keyAndCommentOrder.filter((t) => t !== e)), super.delete(e);
  }
  /**
   * Clear all lines and comments
   */
  clear() {
    this._keyAndCommentOrder.length = 0, this._continuedLine = void 0, this.indent = void 0, this.name = void 0, this.nameKey = void 0, super.clear();
  }
  /**
   * @returns {string} Returns the stanza as a string fit for writing to a ra
   * file. Original leading indent is preserved. It may not be the same as the
   * input stanza as lines that were joined with `\` in the input will be output
   * as a single line and all comments will have the same indentations as the
   * rest of the stanza. Comments between joined lines will move before that
   * line.
   */
  toString() {
    if (this.size === 0)
      return "";
    const e = [];
    return this._keyAndCommentOrder.forEach((t) => {
      t.startsWith("#") ? e.push(`${this.indent}${t}`) : e.push(`${this.indent}${t} ${this.get(t)}`.trimEnd());
    }), `${e.join(`
`)}
`;
  }
}
class f extends Map {
  constructor(e, t = { checkIndent: !0 }) {
    super();
    const { checkIndent: n } = t;
    this._checkIndent = n;
    let i;
    typeof e == "string" ? i = e.trimEnd().split(/(?:[\t ]*\r?\n){2,}/) : e ? i = e : i = [], this._stanzaAndCommentOrder = [], i.forEach((s) => {
      this.add(s);
    });
  }
  /**
   * Add a single stanza to the file
   * @param {string} stanza A single stanza
   * @returns {RaFile} The RaFile object
   */
  add(e) {
    if (e === "")
      throw new Error("Invalid stanza, was empty");
    if (e.trim().startsWith("#")) {
      const n = e.trimEnd().split(/\r?\n/).map((i) => i.trim());
      if (n.every((i) => i.startsWith("#")))
        return this._stanzaAndCommentOrder.push(n.join(`
`)), this;
    }
    const t = new m(e, { checkIndent: this._checkIndent });
    if (!this.nameKey)
      this.nameKey = t.nameKey;
    else if (t.nameKey !== this.nameKey)
      throw new Error(`The first line in each stanza must have the same key. Saw both ${this.nameKey} and ${t.nameKey}`);
    if (!t.name)
      throw new Error(`No stanza name: ${t.name}`);
    if (this.has(t.name))
      throw new Error(`Got duplicate stanza name: ${t.name}`);
    return this._stanzaAndCommentOrder.push(t.name), super.set(t.name, t);
  }
  /**
   * Use `add()` if possible instead of this method. If using this, be aware
   * that no checks are made for comments, empty stanzas, duplicate keys, etc.
   * @param {string} key The key of the RaFile stanza
   * @param {RaStanza} value The RaFile stanza used to replace the prior one
   */
  update(e, t) {
    if (!(t instanceof m))
      throw new Error(`Value of ${e} is not an RaStanza`);
    super.set(e, t);
  }
  /**
   * Delete a stanza
   * @param {string} stanza The name of the stanza to delete (the value in its
   * first key-value pair)
   * @returns {boolean} true if the deleted stanza existed, false if it did not
   */
  delete(e) {
    return this._stanzaAndCommentOrder.includes(e) && (this._stanzaAndCommentOrder = this._stanzaAndCommentOrder.filter((t) => t !== e)), super.delete(e);
  }
  /**
   * Clear all stanzas and comments
   */
  clear() {
    this._stanzaAndCommentOrder.length = 0, this.nameKey = void 0, super.clear();
  }
  /**
   * @returns {string} Returns the stanza as a string fit for writing to a ra
   * file. Original leading indent is preserved. It may not be the same as the
   * input stanza as lines that were joined with `\` in the input will be output
   *  as a single line and all comments will have the same indentations as the
   * rest of the stanza. Comments between joined lines will move before that
   * line.
   */
  toString() {
    if (this.size === 0)
      return "";
    const e = [];
    return this._stanzaAndCommentOrder.forEach((t) => {
      if (t.startsWith("#"))
        e.push(`${t}
`);
      else {
        const n = this.get(t);
        n && e.push(n.toString());
      }
    }), e.join(`
`);
  }
}
class w extends f {
  constructor(e) {
    if (super(e, { checkIndent: !1 }), this.nameKey !== "track")
      throw new Error(`trackDb has "${this.nameKey}" instead of "track" as the first line in each track`);
    this.forEach((t, n) => {
      var i;
      const s = Array.from(t.keys()), r = [];
      if (["track", "shortLabel"].forEach((a) => {
        s.includes(a) || r.push(a);
      }), r.length > 0)
        throw new Error(`Track ${n} is missing required key(s): ${r.join(", ")}`);
      const u = [
        "superTrack",
        "compositeTrack",
        "container",
        "view"
      ];
      if (!s.some((a) => u.includes(a))) {
        if (!s.includes("bigDataUrl"))
          throw new Error(`Track ${n} is missing required key "bigDataUrl"`);
        if (!s.includes("type")) {
          const a = this.settings(n);
          if (!Array.from(a.keys()).includes("type"))
            throw new Error(`Neither track ${n} nor any of its parent tracks have the required key "type"`);
        }
      }
      let c = "", o = n;
      do
        o = (i = this.get(o)) === null || i === void 0 ? void 0 : i.get("parent"), o && ([o] = o.split(" "), c += "    ");
      while (o);
      const l = this.get(n);
      l && (l.indent = c, this.set(n, l));
    });
  }
  /**
   * Gets all track entries including those of parent tracks, with closer
   * entries overriding more distant ones
   * @param {string} trackName The name of a track
   * @throws {Error} Throws if track name does not exist in the trackDb
   */
  settings(e) {
    var t;
    if (!this.has(e))
      throw new Error(`Track ${e} does not exist`);
    const n = [e];
    let i = e;
    do
      i = (t = this.get(i)) === null || t === void 0 ? void 0 : t.get("parent"), i && n.push(i);
    while (i);
    const s = /* @__PURE__ */ new Map();
    return n.reverse(), n.forEach((r) => {
      var h;
      (h = this.get(r)) === null || h === void 0 || h.forEach((u, c) => {
        s.set(c, u);
      });
    }), s;
  }
}
class g extends m {
  constructor(e) {
    if (super(e), this.nameKey !== "hub")
      throw new Error('Hub file must begin with a line like "hub <hub_name>"');
    const t = [
      "hub",
      "shortLabel",
      "longLabel",
      "genomesFile",
      "email",
      "descriptionUrl"
    ], n = [];
    if (this.forEach((s, r) => {
      t.includes(r) || n.push(r);
    }), n.length > 0)
      throw new Error(`Hub file has invalid entr${n.length === 1 ? "y" : "ies"}: ${n.join(", ")}`);
    const i = [];
    if (t.forEach((s) => {
      s !== "descriptionUrl" && !this.get(s) && i.push(s);
    }), i.length > 0)
      throw new Error(`Hub file is missing required entr${i.length === 1 ? "y" : "ies"}: ${i.join(", ")}`);
  }
}
class y extends f {
  constructor(e) {
    if (super(e), this.nameKey !== "genome")
      throw new Error('Genomes file must begin with a line like "genome <genome_name>"');
    const t = [
      "genome",
      "trackDb"
      // 'twoBitPath',
      // 'groups',
    ];
    this.forEach((n, i) => {
      const s = [];
      if (t.forEach((r) => {
        n.get(r) || s.push(r);
      }), s.length > 0)
        throw new Error(`Genomes file entry ${i} is missing required entr${s.length === 1 ? "y" : "ies"}: ${s.join(", ")}`);
    });
  }
}
export {
  y as GenomesFile,
  g as HubFile,
  f as RaFile,
  m as RaStanza,
  w as TrackDbFile
};
