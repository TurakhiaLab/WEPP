import { aJ as v } from "./index-8ee3baa4.js";
import { r as R, s as D } from "./index-aa7a8a76.js";
var B = R();
function p(t) {
  return t.replace(/%([0-9A-Fa-f]{2})/g, (e, s) => String.fromCharCode(parseInt(s, 16)));
}
function F(t, e) {
  return String(e).replace(t, (s) => `%${s.charCodeAt(0).toString(16).toUpperCase().padStart(2, "0")}`);
}
function d(t) {
  return F(/[\n;\r\t=%&,\x00-\x1f\x7f-\xff]/g, t);
}
function u(t) {
  return F(/[\n\r\t%\x00-\x1f\x7f-\xff]/g, t);
}
function k(t) {
  if (!(t && t.length) || t === ".")
    return {};
  const e = {};
  return t.replace(/\r?\n$/, "").split(";").forEach((s) => {
    const r = s.split("=", 2);
    if (!(r[1] && r[1].length))
      return;
    r[0] = r[0].trim();
    let n = e[r[0].trim()];
    n || (n = [], e[r[0]] = n), n.push(...r[1].split(",").map((i) => i.trim()).map(p));
  }), e;
}
function A(t) {
  const e = t.split("	").map((r) => r === "." || r === "" ? null : r);
  return {
    seq_id: e[0] && p(e[0]),
    source: e[1] && p(e[1]),
    type: e[2] && p(e[2]),
    start: e[3] === null ? null : parseInt(e[3], 10),
    end: e[4] === null ? null : parseInt(e[4], 10),
    score: e[5] === null ? null : parseFloat(e[5]),
    strand: e[6],
    phase: e[7],
    attributes: e[8] === null ? null : k(e[8])
  };
}
function x(t) {
  const e = /^\s*##\s*(\S+)\s*(.*)/.exec(t);
  if (!e)
    return null;
  const [, s] = e;
  let [, , r] = e;
  const n = { directive: s };
  if (r.length && (r = r.replace(/\r?\n$/, ""), n.value = r), s === "sequence-region") {
    const i = r.split(/\s+/, 3);
    return {
      ...n,
      seq_id: i[0],
      start: i[1] && i[1].replace(/\D/g, ""),
      end: i[2] && i[2].replace(/\D/g, "")
    };
  } else if (s === "genome-build") {
    const [i, a] = r.split(/\s+/, 2);
    return {
      ...n,
      source: i,
      buildName: a
    };
  }
  return n;
}
function q(t) {
  const e = [];
  return Object.entries(t).forEach(([s, r]) => {
    if (!r)
      return;
    let n;
    r.hasOwnProperty("toString") ? n = d(r.toString()) : Array.isArray(r) ? n = r.map(d).join(",") : n = d(r), e.push(`${d(s)}=${n}`);
  }), e.length ? e.join(";") : ".";
}
function O(t, e) {
  const s = t.attributes === null || t.attributes === void 0 ? "." : q(t.attributes), n = `${[
    t.seq_id === null ? "." : u(t.seq_id),
    t.source === null ? "." : u(t.source),
    t.type === null ? "." : u(t.type),
    t.start === null ? "." : u(t.start),
    t.end === null ? "." : u(t.end),
    t.score === null ? "." : u(t.score),
    t.strand === null ? "." : u(t.strand),
    t.phase === null ? "." : u(t.phase),
    s
  ].join("	")}
`;
  return e[n] ? "" : (e[n] = !0, n);
}
function _(t, e) {
  if (Array.isArray(t))
    return t.map((r) => _(r, e)).join("");
  const s = [O(t, e)];
  return j(t) && s.push(...t.child_features.map((r) => _(r, e)), ...t.derived_features.map((r) => _(r, e))), s.join("");
}
function L(t) {
  return _(t, {});
}
function E(t) {
  let e = `##${t.directive}`;
  return t.value && (e += ` ${t.value}`), e += `
`, e;
}
function P(t) {
  return `# ${t.comment}
`;
}
function S(t) {
  return `>${t.id}${t.description ? ` ${t.description}` : ""}
${t.sequence}
`;
}
function b(t) {
  function e(s) {
    return "attributes" in s ? L(s) : "directive" in s ? E(s) : "sequence" in s ? S(s) : "comment" in s ? P(s) : `# (invalid item found during format)
`;
  }
  return Array.isArray(t) ? t.map(e) : e(t);
}
function j(t) {
  return t.child_features !== void 0 && t.derived_features !== void 0;
}
const I = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  escape: d,
  escapeColumn: u,
  formatAttributes: q,
  formatComment: P,
  formatDirective: E,
  formatFeature: L,
  formatItem: b,
  formatSequence: S,
  parseAttributes: k,
  parseDirective: x,
  parseFeature: A,
  unescape: p
}, Symbol.toStringTag, { value: "Module" })), g = {
  Parent: "child_features",
  Derives_from: "derived_features"
};
class y {
  constructor(e) {
    this.seqCallback = e, this.currentSequence = void 0;
  }
  addLine(e) {
    const s = /^>\s*(\S+)\s*(.*)/.exec(e);
    s ? (this._flush(), this.currentSequence = { id: s[1], sequence: "" }, s[2] && (this.currentSequence.description = s[2].trim())) : this.currentSequence && /\S/.test(e) && (this.currentSequence.sequence += e.replace(/\s/g, ""));
  }
  _flush() {
    this.currentSequence && this.seqCallback(this.currentSequence);
  }
  finish() {
    this._flush();
  }
}
class $ {
  constructor(e) {
    this.fastaParser = void 0, this.eof = !1, this.lineNumber = 0, this._underConstructionTopLevel = [], this._underConstructionById = {}, this._completedReferences = {}, this._underConstructionOrphans = {};
    const s = () => {
    };
    this.featureCallback = e.featureCallback || s, this.endCallback = e.endCallback || s, this.commentCallback = e.commentCallback || s, this.errorCallback = e.errorCallback || s, this.directiveCallback = e.directiveCallback || s, this.sequenceCallback = e.sequenceCallback || s, this.disableDerivesFromReferences = e.disableDerivesFromReferences || !1, this.bufferSize = e.bufferSize === void 0 ? 1e3 : e.bufferSize;
  }
  addLine(e) {
    if (this.fastaParser) {
      this.fastaParser.addLine(e);
      return;
    }
    if (this.eof)
      return;
    if (this.lineNumber += 1, /^\s*[^#\s>]/.test(e)) {
      this._bufferLine(e);
      return;
    }
    const s = /^\s*(#+)(.*)/.exec(e);
    if (s) {
      const [, r] = s;
      let [, , n] = s;
      if (r.length === 3)
        this._emitAllUnderConstructionFeatures();
      else if (r.length === 2) {
        const i = x(e);
        i && (i.directive === "FASTA" ? (this._emitAllUnderConstructionFeatures(), this.eof = !0, this.fastaParser = new y(this.sequenceCallback)) : this._emitItem(i));
      } else
        n = n.replace(/\s*/, ""), this._emitItem({ comment: n });
    } else if (!/^\s*$/.test(e))
      if (/^\s*>/.test(e))
        this._emitAllUnderConstructionFeatures(), this.eof = !0, this.fastaParser = new y(this.sequenceCallback), this.fastaParser.addLine(e);
      else {
        const r = e.replace(/\r?\n?$/g, "");
        throw new Error(`GFF3 parse error.  Cannot parse '${r}'.`);
      }
  }
  finish() {
    this._emitAllUnderConstructionFeatures(), this.fastaParser && this.fastaParser.finish(), this.endCallback();
  }
  _emitItem(e) {
    Array.isArray(e) ? this.featureCallback(e) : "directive" in e ? this.directiveCallback(e) : "comment" in e && this.commentCallback(e);
  }
  _enforceBufferSizeLimit(e = 0) {
    const s = (r) => {
      r && Array.isArray(r) && r[0].attributes && r[0].attributes.ID && r[0].attributes.ID[0] && (r[0].attributes.ID.forEach((i) => {
        delete this._underConstructionById[i], delete this._completedReferences[i];
      }), r.forEach((i) => {
        i.child_features && i.child_features.forEach((a) => s(a)), i.derived_features && i.derived_features.forEach((a) => s(a));
      }));
    };
    for (; this._underConstructionTopLevel.length + e > this.bufferSize; ) {
      const r = this._underConstructionTopLevel.shift();
      r && (this._emitItem(r), s(r));
    }
  }
  /**
   * return all under-construction features, called when we know
   * there will be no additional data to attach to them
   */
  _emitAllUnderConstructionFeatures() {
    if (this._underConstructionTopLevel.forEach(this._emitItem.bind(this)), this._underConstructionTopLevel = [], this._underConstructionById = {}, this._completedReferences = {}, Array.from(Object.values(this._underConstructionOrphans)).length)
      throw new Error(`some features reference other features that do not exist in the file (or in the same '###' scope). ${Object.keys(this._underConstructionOrphans)}`);
  }
  // do the right thing with a newly-parsed feature line
  _bufferLine(e) {
    var s, r, n;
    const a = {
      ...A(e),
      child_features: [],
      derived_features: []
    }, o = ((s = a.attributes) === null || s === void 0 ? void 0 : s.ID) || [], c = ((r = a.attributes) === null || r === void 0 ? void 0 : r.Parent) || [], h = this.disableDerivesFromReferences ? [] : ((n = a.attributes) === null || n === void 0 ? void 0 : n.Derives_from) || [];
    if (!o.length && !c.length && !h.length) {
      this._emitItem([a]);
      return;
    }
    let l;
    o.forEach((m) => {
      const f = this._underConstructionById[m];
      f ? (f[f.length - 1].type !== a.type && this._parseError(`multi-line feature "${m}" has inconsistent types: "${a.type}", "${f[f.length - 1].type}"`), f.push(a), l = f) : (l = [a], this._enforceBufferSizeLimit(1), !c.length && !h.length && this._underConstructionTopLevel.push(l), this._underConstructionById[m] = l, this._resolveReferencesTo(l, m));
    }), this._resolveReferencesFrom(l || [a], { Parent: c, Derives_from: h }, o);
  }
  _resolveReferencesTo(e, s) {
    const r = this._underConstructionOrphans[s];
    r && (e.forEach((n) => {
      n.child_features.push(...r.Parent);
    }), e.forEach((n) => {
      n.derived_features.push(...r.Derives_from);
    }), delete this._underConstructionOrphans[s]);
  }
  _parseError(e) {
    this.eof = !0, this.errorCallback(`${this.lineNumber}: ${e}`);
  }
  _resolveReferencesFrom(e, s, r) {
    function n(i, a, o) {
      let c = i[a];
      c || (c = {}, i[a] = c);
      const h = c[o] || !1;
      return c[o] = !0, h;
    }
    s.Parent.forEach((i) => {
      const a = this._underConstructionById[i];
      if (a) {
        const o = g.Parent;
        r.filter((c) => n(this._completedReferences, c, `Parent,${i}`)).length || a.forEach((c) => {
          c[o].push(e);
        });
      } else {
        let o = this._underConstructionOrphans[i];
        o || (o = {
          Parent: [],
          Derives_from: []
        }, this._underConstructionOrphans[i] = o), o.Parent.push(e);
      }
    }), s.Derives_from.forEach((i) => {
      const a = this._underConstructionById[i];
      if (a) {
        const o = g.Derives_from;
        r.filter((c) => n(this._completedReferences, c, `Derives_from,${i}`)).length || a.forEach((c) => {
          c[o].push(e);
        });
      } else {
        let o = this._underConstructionOrphans[i];
        o || (o = {
          Parent: [],
          Derives_from: []
        }, this._underConstructionOrphans[i] = o), o.Derives_from.push(e);
      }
    });
  }
}
function C(t) {
  v && v.nextTick ? v.nextTick(t) : t();
}
function w(t) {
  const e = {
    encoding: "utf8",
    parseFeatures: !0,
    parseDirectives: !1,
    parseSequences: !0,
    parseComments: !1,
    bufferSize: 1e3,
    disableDerivesFromReferences: !1,
    ...t
  };
  return t.parseAll && (e.parseFeatures = !0, e.parseDirectives = !0, e.parseComments = !0, e.parseSequences = !0), e;
}
class M extends D.Transform {
  constructor(e = {}) {
    super({ objectMode: !0 }), this.textBuffer = "";
    const s = w(e);
    this.encoding = e.encoding || "utf8", this.decoder = new B.StringDecoder();
    const r = this.push.bind(this);
    this.parser = new $({
      featureCallback: s.parseFeatures ? r : void 0,
      directiveCallback: s.parseDirectives ? r : void 0,
      commentCallback: s.parseComments ? r : void 0,
      sequenceCallback: s.parseSequences ? r : void 0,
      errorCallback: (n) => this.emit("error", n),
      bufferSize: s.bufferSize,
      disableDerivesFromReferences: s.disableDerivesFromReferences
    });
  }
  _addLine(e) {
    e && this.parser.addLine(e);
  }
  _nextText(e) {
    const s = (this.textBuffer + e).split(/\r?\n/);
    this.textBuffer = s.pop() || "", s.forEach((r) => this._addLine(r));
  }
  _transform(e, s, r) {
    this._nextText(this.decoder.write(e)), C(r);
  }
  _flush(e) {
    this.decoder.end && this._nextText(this.decoder.end()), this.textBuffer != null && this._addLine(this.textBuffer), this.parser.finish(), C(e);
  }
}
function z(t = {}) {
  return new M(t);
}
function U(t, e = {}) {
  if (!t)
    return [];
  const s = w(e), r = [], n = r.push.bind(r), i = new $({
    featureCallback: s.parseFeatures ? n : void 0,
    directiveCallback: s.parseDirectives ? n : void 0,
    commentCallback: s.parseComments ? n : void 0,
    sequenceCallback: s.parseSequences ? n : void 0,
    disableDerivesFromReferences: s.disableDerivesFromReferences || !1,
    bufferSize: 1 / 0,
    errorCallback: (a) => {
      throw a;
    }
  });
  return t.split(/\r?\n/).forEach(i.addLine.bind(i)), i.finish(), r;
}
function V(t) {
  const e = [], s = [];
  t.forEach((n) => {
    "sequence" in n ? s.push(n) : e.push(n);
  });
  let r = e.map(b).join("");
  return s.length && (r += `##FASTA
`, r += s.map(S).join("")), r;
}
class T extends D.Transform {
  constructor(e = {}) {
    super(Object.assign(e, { objectMode: !0 })), this.linesSinceLastSyncMark = 0, this.haveWeEmittedData = !1, this.fastaMode = !1, this.minLinesBetweenSyncMarks = e.minSyncLines || 100, this.insertVersionDirective = e.insertVersionDirective || !1;
  }
  _transform(e, s, r) {
    let n;
    if (!this.haveWeEmittedData && this.insertVersionDirective) {
      const i = Array.isArray(e) ? e[0] : e;
      "directive" in i && i.directive !== "gff-version" && this.push(`##gff-version 3
`);
    }
    if ("sequence" in e && !this.fastaMode && (this.push(`##FASTA
`), this.fastaMode = !0), Array.isArray(e) ? n = e.map(b).join("") : n = b(e), this.push(n), this.linesSinceLastSyncMark >= this.minLinesBetweenSyncMarks)
      this.push(`###
`), this.linesSinceLastSyncMark = 0;
    else {
      let i = 0;
      for (let a = 0; a < n.length; a += 1)
        n[a] === `
` && (i += 1);
      this.linesSinceLastSyncMark += i;
    }
    this.haveWeEmittedData = !0, C(r);
  }
}
function N(t = {}) {
  return new T(t);
}
function W(t, e, s = {}) {
  const r = {
    insertVersionDirective: !0,
    ...s
  };
  return new Promise((n, i) => {
    t.pipe(new T(r)).on("end", () => n(null)).on("error", i).pipe(e);
  });
}
const H = {
  parseStream: z,
  parseStringSync: U,
  formatSync: V,
  formatStream: N,
  formatFile: W,
  util: I
};
export {
  H as g
};
