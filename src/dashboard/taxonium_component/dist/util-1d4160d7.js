import { Q as xt, t as Ie } from "./index-8ee3baa4.js";
function Nt(s, m) {
  function t() {
    this.constructor = s;
  }
  t.prototype = m.prototype, s.prototype = new t();
}
function j(s, m, t, d) {
  this.message = s, this.expected = m, this.found = t, this.location = d, this.name = "SyntaxError", typeof Error.captureStackTrace == "function" && Error.captureStackTrace(this, j);
}
Nt(j, Error);
j.buildMessage = function(s, m) {
  var t = {
    literal: function(a) {
      return '"' + w(a.text) + '"';
    },
    class: function(a) {
      var f = "", p;
      for (p = 0; p < a.parts.length; p++)
        f += a.parts[p] instanceof Array ? E(a.parts[p][0]) + "-" + E(a.parts[p][1]) : E(a.parts[p]);
      return "[" + (a.inverted ? "^" : "") + f + "]";
    },
    any: function(a) {
      return "any character";
    },
    end: function(a) {
      return "end of input";
    },
    other: function(a) {
      return a.description;
    }
  };
  function d(a) {
    return a.charCodeAt(0).toString(16).toUpperCase();
  }
  function w(a) {
    return a.replace(/\\/g, "\\\\").replace(/"/g, '\\"').replace(/\0/g, "\\0").replace(/\t/g, "\\t").replace(/\n/g, "\\n").replace(/\r/g, "\\r").replace(/[\x00-\x0F]/g, function(f) {
      return "\\x0" + d(f);
    }).replace(/[\x10-\x1F\x7F-\x9F]/g, function(f) {
      return "\\x" + d(f);
    });
  }
  function E(a) {
    return a.replace(/\\/g, "\\\\").replace(/\]/g, "\\]").replace(/\^/g, "\\^").replace(/-/g, "\\-").replace(/\0/g, "\\0").replace(/\t/g, "\\t").replace(/\n/g, "\\n").replace(/\r/g, "\\r").replace(/[\x00-\x0F]/g, function(f) {
      return "\\x0" + d(f);
    }).replace(/[\x10-\x1F\x7F-\x9F]/g, function(f) {
      return "\\x" + d(f);
    });
  }
  function h(a) {
    return t[a.type](a);
  }
  function k(a) {
    var f = new Array(a.length), p, v;
    for (p = 0; p < a.length; p++)
      f[p] = h(a[p]);
    if (f.sort(), f.length > 0) {
      for (p = 1, v = 1; p < f.length; p++)
        f[p - 1] !== f[p] && (f[v] = f[p], v++);
      f.length = v;
    }
    switch (f.length) {
      case 1:
        return f[0];
      case 2:
        return f[0] + " or " + f[1];
      default:
        return f.slice(0, -1).join(", ") + ", or " + f[f.length - 1];
    }
  }
  function y(a) {
    return a ? '"' + w(a) + '"' : "end of input";
  }
  return "Expected " + k(s) + " but " + y(m) + " found.";
};
function qt(s, m) {
  m = m !== void 0 ? m : {};
  var t = {}, d = { declaration: Te }, w = Te, E = "(", h = b("(", !1), k = ")", y = b(")", !1), a = function(e, o, i, n) {
    return { type: e, name: o, comment: i, fields: n };
  }, f = "simple", p = b("simple", !1), v = "object", A = b("object", !1), q = "table", X = b("table", !1), I = "auto", O = b("auto", !1), G = "primary", Y = b("primary", !1), L = "index", H = b("index", !1), z = "unique", V = b("unique", !1), U = function(e, o) {
    return o;
  }, Ue = function(e, o) {
    return e.name && o.unshift(e), o;
  }, je = "#", Me = b("#", !1), J = ";", K = b(";", !1), Oe = function(e, o, i) {
    return { type: e, name: o, comment: i };
  }, Le = "[", Ve = b("[", !1), Ze = "]", Qe = b("]", !1), We = function(e, o, i, n) {
    return { type: e, size: o, name: i, comment: n };
  }, Xe = function(e, o, i, n) {
    return { type: e, vals: o, name: i, comment: n };
  }, ie = ",", ne = b(",", !1), Ye = function(e, o) {
    return o.unshift(e), o;
  }, ae = "int", He = b("int", !1), ce = "uint", Je = b("uint", !1), fe = "short", Ke = b("short", !1), le = "ushort", et = b("ushort", !1), ue = "byte", tt = b("byte", !1), he = "ubyte", rt = b("ubyte", !1), me = "float", ot = b("float", !1), pe = "char", st = b("char", !1), ge = "string", it = b("string", !1), de = "lstring", nt = b("lstring", !1), be = "enum", at = b("enum", !1), Se = "double", ct = b("double", !1), $e = "bigint", ft = b("bigint", !1), ye = "set", lt = b("set", !1), ut = function(e, o) {
    return e + " " + o;
  }, ht = /^[a-zA-Z_]/, mt = M([["a", "z"], ["A", "Z"], "_"], !1, !1), Ee = /^[a-zA-Z0-9_]/, ke = M([["a", "z"], ["A", "Z"], ["0", "9"], "_"], !1, !1), pt = function(e) {
    return qe();
  }, Ce = /^[^\n\r]/, we = M([`
`, "\r"], !0, !1), gt = function(e) {
    return e.join("").replace(/^"/, "").replace(/"$/, "");
  }, dt = De("integer"), ve = /^[0-9]/, Ae = M([["0", "9"]], !1, !1), bt = function() {
    return parseInt(qe(), 10);
  }, St = De("whitespace"), xe = /^[ \t\n\r]/, Ne = M([" ", "	", `
`, "\r"], !1, !1), r = 0, C = 0, Z = [{ line: 1, column: 1 }], D = 0, ee = [], l = 0, Q;
  if ("startRule" in m) {
    if (!(m.startRule in d))
      throw new Error(`Can't start parsing from rule "` + m.startRule + '".');
    w = d[m.startRule];
  }
  function qe() {
    return s.substring(C, r);
  }
  function b(e, o) {
    return { type: "literal", text: e, ignoreCase: o };
  }
  function M(e, o, i) {
    return { type: "class", parts: e, inverted: o, ignoreCase: i };
  }
  function $t() {
    return { type: "end" };
  }
  function De(e) {
    return { type: "other", description: e };
  }
  function Fe(e) {
    var o = Z[e], i;
    if (o)
      return o;
    for (i = e - 1; !Z[i]; )
      i--;
    for (o = Z[i], o = {
      line: o.line,
      column: o.column
    }; i < e; )
      s.charCodeAt(i) === 10 ? (o.line++, o.column = 1) : o.column++, i++;
    return Z[e] = o, o;
  }
  function Re(e, o) {
    var i = Fe(e), n = Fe(o);
    return {
      start: {
        offset: e,
        line: i.line,
        column: i.column
      },
      end: {
        offset: o,
        line: n.line,
        column: n.column
      }
    };
  }
  function u(e) {
    r < D || (r > D && (D = r, ee = []), ee.push(e));
  }
  function yt(e, o, i) {
    return new j(j.buildMessage(e, o), e, o, i);
  }
  function Te() {
    var e, o, i, n, c, S, $, x, B, F, _, R, P, T;
    return e = r, o = g(), o !== t ? (i = ze(), i !== t ? (n = g(), n !== t ? (c = Be(), c !== t ? (S = g(), S !== t ? ($ = W(), $ !== t ? (x = g(), x !== t ? (s.charCodeAt(r) === 40 ? (B = E, r++) : (B = t, l === 0 && u(h)), B !== t ? (F = g(), F !== t ? (_ = Et(), _ !== t ? (R = g(), R !== t ? (s.charCodeAt(r) === 41 ? (P = k, r++) : (P = t, l === 0 && u(y)), P !== t ? (T = g(), T !== t ? (C = e, o = a(i, c, $, _), e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t), e;
  }
  function ze() {
    var e;
    return s.substr(r, 6) === f ? (e = f, r += 6) : (e = t, l === 0 && u(p)), e === t && (s.substr(r, 6) === v ? (e = v, r += 6) : (e = t, l === 0 && u(A)), e === t && (s.substr(r, 5) === q ? (e = q, r += 5) : (e = t, l === 0 && u(X)))), e;
  }
  function Be() {
    var e, o, i, n;
    return e = N(), e === t && (e = r, o = N(), o !== t ? (i = _e(), i !== t ? (o = [o, i], e = o) : (r = e, e = t)) : (r = e, e = t), e === t && (e = r, o = N(), o !== t ? (s.substr(r, 4) === I ? (i = I, r += 4) : (i = t, l === 0 && u(O)), i !== t ? (o = [o, i], e = o) : (r = e, e = t)) : (r = e, e = t), e === t && (e = r, o = N(), o !== t ? (i = _e(), i !== t ? (s.substr(r, 4) === I ? (n = I, r += 4) : (n = t, l === 0 && u(O)), n !== t ? (o = [o, i, n], e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)))), e;
  }
  function _e() {
    var e;
    return s.substr(r, 7) === G ? (e = G, r += 7) : (e = t, l === 0 && u(Y)), e === t && (s.substr(r, 5) === L ? (e = L, r += 5) : (e = t, l === 0 && u(H)), e === t && (s.substr(r, 6) === z ? (e = z, r += 6) : (e = t, l === 0 && u(V)))), e;
  }
  function W() {
    var e;
    return e = Pe(), e === t && (e = g()), e;
  }
  function Et() {
    var e, o, i, n, c, S, $;
    if (e = r, o = te(), o !== t)
      if (i = g(), i !== t) {
        for (n = [], c = r, S = g(), S !== t ? ($ = te(), $ !== t ? (C = c, S = U(o, $), c = S) : (r = c, c = t)) : (r = c, c = t); c !== t; )
          n.push(c), c = r, S = g(), S !== t ? ($ = te(), $ !== t ? (C = c, S = U(o, $), c = S) : (r = c, c = t)) : (r = c, c = t);
        n !== t ? (c = g(), c !== t ? (C = e, o = Ue(o, n), e = o) : (r = e, e = t)) : (r = e, e = t);
      } else
        r = e, e = t;
    else
      r = e, e = t;
    return e;
  }
  function kt() {
    var e;
    return s.charCodeAt(r) === 35 ? (e = je, r++) : (e = t, l === 0 && u(Me)), e;
  }
  function Ct() {
    var e, o, i, n, c;
    return e = r, o = g(), o !== t ? (i = kt(), i !== t ? (n = Pe(), n !== t ? (c = g(), c !== t ? (o = [o, i, n, c], e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t), e;
  }
  function te() {
    var e, o, i, n, c, S, $, x, B, F, _, R, P, T;
    return e = r, o = re(), o !== t ? (i = g(), i !== t ? (n = N(), n !== t ? (c = g(), c !== t ? (s.charCodeAt(r) === 59 ? (S = J, r++) : (S = t, l === 0 && u(K)), S !== t ? ($ = g(), $ !== t ? (x = W(), x !== t ? (C = e, o = Oe(o, n, x), e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t), e === t && (e = r, o = re(), o !== t ? (i = g(), i !== t ? (s.charCodeAt(r) === 91 ? (n = Le, r++) : (n = t, l === 0 && u(Ve)), n !== t ? (c = g(), c !== t ? (S = vt(), S !== t ? ($ = g(), $ !== t ? (s.charCodeAt(r) === 93 ? (x = Ze, r++) : (x = t, l === 0 && u(Qe)), x !== t ? (B = g(), B !== t ? (F = N(), F !== t ? (_ = g(), _ !== t ? (s.charCodeAt(r) === 59 ? (R = J, r++) : (R = t, l === 0 && u(K)), R !== t ? (P = g(), P !== t ? (T = W(), T !== t ? (C = e, o = We(o, S, F, T), e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t), e === t && (e = r, o = re(), o !== t ? (i = g(), i !== t ? (s.charCodeAt(r) === 40 ? (n = E, r++) : (n = t, l === 0 && u(h)), n !== t ? (c = g(), c !== t ? (S = wt(), S !== t ? ($ = g(), $ !== t ? (s.charCodeAt(r) === 41 ? (x = k, r++) : (x = t, l === 0 && u(y)), x !== t ? (B = g(), B !== t ? (F = N(), F !== t ? (_ = g(), _ !== t ? (s.charCodeAt(r) === 59 ? (R = J, r++) : (R = t, l === 0 && u(K)), R !== t ? (P = g(), P !== t ? (T = W(), T !== t ? (C = e, o = Xe(o, S, F, T), e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t), e === t && (e = Ct()))), e;
  }
  function wt() {
    var e, o, i, n, c, S, $;
    if (e = r, o = N(), o !== t) {
      for (i = [], n = r, s.charCodeAt(r) === 44 ? (c = ie, r++) : (c = t, l === 0 && u(ne)), c !== t ? (S = g(), S !== t ? ($ = N(), $ !== t ? (C = n, c = U(o, $), n = c) : (r = n, n = t)) : (r = n, n = t)) : (r = n, n = t); n !== t; )
        i.push(n), n = r, s.charCodeAt(r) === 44 ? (c = ie, r++) : (c = t, l === 0 && u(ne)), c !== t ? (S = g(), S !== t ? ($ = N(), $ !== t ? (C = n, c = U(o, $), n = c) : (r = n, n = t)) : (r = n, n = t)) : (r = n, n = t);
      i !== t ? (C = e, o = Ye(o, i), e = o) : (r = e, e = t);
    } else
      r = e, e = t;
    return e;
  }
  function re() {
    var e, o, i, n;
    return s.substr(r, 3) === ae ? (e = ae, r += 3) : (e = t, l === 0 && u(He)), e === t && (s.substr(r, 4) === ce ? (e = ce, r += 4) : (e = t, l === 0 && u(Je)), e === t && (s.substr(r, 5) === fe ? (e = fe, r += 5) : (e = t, l === 0 && u(Ke)), e === t && (s.substr(r, 6) === le ? (e = le, r += 6) : (e = t, l === 0 && u(et)), e === t && (s.substr(r, 4) === ue ? (e = ue, r += 4) : (e = t, l === 0 && u(tt)), e === t && (s.substr(r, 5) === he ? (e = he, r += 5) : (e = t, l === 0 && u(rt)), e === t && (s.substr(r, 5) === me ? (e = me, r += 5) : (e = t, l === 0 && u(ot)), e === t && (s.substr(r, 4) === pe ? (e = pe, r += 4) : (e = t, l === 0 && u(st)), e === t && (s.substr(r, 6) === ge ? (e = ge, r += 6) : (e = t, l === 0 && u(it)), e === t && (s.substr(r, 7) === de ? (e = de, r += 7) : (e = t, l === 0 && u(nt)), e === t && (s.substr(r, 4) === be ? (e = be, r += 4) : (e = t, l === 0 && u(at)), e === t && (s.substr(r, 6) === Se ? (e = Se, r += 6) : (e = t, l === 0 && u(ct)), e === t && (s.substr(r, 6) === $e ? (e = $e, r += 6) : (e = t, l === 0 && u(ft)), e === t && (s.substr(r, 3) === ye ? (e = ye, r += 3) : (e = t, l === 0 && u(lt)), e === t && (e = r, o = ze(), o !== t ? (i = g(), i !== t ? (n = Be(), n !== t ? (C = e, o = ut(o, n), e = o) : (r = e, e = t)) : (r = e, e = t)) : (r = e, e = t))))))))))))))), e;
  }
  function vt() {
    var e;
    return e = At(), e === t && (e = N()), e;
  }
  function N() {
    var e, o, i, n, c;
    if (e = r, o = r, ht.test(s.charAt(r)) ? (i = s.charAt(r), r++) : (i = t, l === 0 && u(mt)), i !== t) {
      for (n = [], Ee.test(s.charAt(r)) ? (c = s.charAt(r), r++) : (c = t, l === 0 && u(ke)); c !== t; )
        n.push(c), Ee.test(s.charAt(r)) ? (c = s.charAt(r), r++) : (c = t, l === 0 && u(ke));
      n !== t ? (i = [i, n], o = i) : (r = o, o = t);
    } else
      r = o, o = t;
    return o !== t && (C = e, o = pt()), e = o, e;
  }
  function Pe() {
    var e, o, i;
    for (e = r, o = [], Ce.test(s.charAt(r)) ? (i = s.charAt(r), r++) : (i = t, l === 0 && u(we)); i !== t; )
      o.push(i), Ce.test(s.charAt(r)) ? (i = s.charAt(r), r++) : (i = t, l === 0 && u(we));
    return o !== t && (C = e, o = gt(o)), e = o, e;
  }
  function At() {
    var e, o, i, n;
    if (l++, e = r, o = g(), o !== t) {
      if (i = [], ve.test(s.charAt(r)) ? (n = s.charAt(r), r++) : (n = t, l === 0 && u(Ae)), n !== t)
        for (; n !== t; )
          i.push(n), ve.test(s.charAt(r)) ? (n = s.charAt(r), r++) : (n = t, l === 0 && u(Ae));
      else
        i = t;
      i !== t ? (C = e, o = bt(), e = o) : (r = e, e = t);
    } else
      r = e, e = t;
    return l--, e === t && (o = t, l === 0 && u(dt)), e;
  }
  function g() {
    var e, o;
    for (l++, e = [], xe.test(s.charAt(r)) ? (o = s.charAt(r), r++) : (o = t, l === 0 && u(Ne)); o !== t; )
      e.push(o), xe.test(s.charAt(r)) ? (o = s.charAt(r), r++) : (o = t, l === 0 && u(Ne));
    return l--, e === t && (o = t, l === 0 && u(St)), e;
  }
  if (Q = w(), Q !== t && r === s.length)
    return Q;
  throw Q !== t && r < s.length && u($t()), yt(ee, D < s.length ? s.charAt(D) : null, D < s.length ? Re(D, D + 1) : Re(D, D));
}
var Ge = {
  SyntaxError: j,
  parse: qt
};
const Dt = /* @__PURE__ */ xt(Ge), Ft = `table bigChain
"bigChain pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint tSize;         "size of target sequence"
    string qName;       "name of query sequence"
    uint qSize;         "size of query sequence"
    uint qStart;        "start of alignment on query sequence"
    uint qEnd;          "end of alignment on query sequence"
    uint chainScore;    "score from chain"
    )`, Rt = `table bigGenePred
"bigGenePred gene models"
   (
   string chrom;       "Reference sequence chromosome or scaffold"
   uint   chromStart;  "Start position in chromosome"
   uint   chromEnd;    "End position in chromosome"
   string name;        "Name or ID of item, ideally both human readable and unique"
   uint score;         "Score (0-1000)"
   char[1] strand;     "+ or - for strand"
   uint thickStart;    "Start of where display should be thick (start codon)"
   uint thickEnd;      "End of where display should be thick (stop codon)"
   uint reserved;       "RGB value (use R,G,B string in input file)"
   int blockCount;     "Number of blocks"
   int[blockCount] blockSizes; "Comma separated list of block sizes"
   int[blockCount] chromStarts; "Start positions relative to chromStart"
   string name2;       "Alternative/human readable name"
   string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
   string cdsEndStat;   "Status of CDS end annotation (none, unknown, incomplete, or complete)"
   int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
   string type;        "Transcript type"
   string geneName;    "Primary identifier for gene"
   string geneName2;   "Alternative/human readable gene name"
   string geneType;    "Gene type"
   )`, Tt = `table interact
"interaction between two regions"
    (
    string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
    uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region"
    uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
    uint score;          "Score (0-1000)"
    double value;        "Strength of interaction or other data value. Typically basis for score"
    string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
    string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;    "Start position in chromosome of source/lower/this region"
    uint sourceEnd;      "End position in chromosome of source/lower/this region"
    string sourceName;   "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;    "Start position in chromosome of target/upper/this region"
    uint targetEnd;      "End position in chromosome of target/upper/this region"
    string targetName;   "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"

    )`, zt = `table bigLink
"bigLink pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint qStart;        "start of alignment on query sequence"
    )`, Bt = `table bedMaf
"Bed3 with MAF block"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    lstring mafBlock;   "MAF block"
    )`, _t = `table bigNarrowPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    string name;	 "Name given to a region (preferably unique). Use . if no name is assigned"
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
    int   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)`, Pt = `table bigPsl
"bigPsl pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - indicates whether the query aligns to the + or - strand on the reference"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"

    uint    oChromStart;"Start position in other chromosome"
    uint    oChromEnd;  "End position in other chromosome"
    char[1] oStrand;    "+ or -, - means that psl was reversed into BED-compatible coordinates"
    uint    oChromSize; "Size of other chromosome."
    int[blockCount] oChromStarts; "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"

    lstring  oSequence;  "Sequence on other chrom (or empty)"
    string   oCDS;       "CDS in NCBI format"

    uint    chromSize;"Size of target chromosome"

    uint match;        "Number of bases matched."
    uint misMatch; " Number of bases that don't match "
    uint repMatch; " Number of bases that match but are part of repeats "
    uint nCount;   " Number of 'N' bases "
    uint seqType;    "0=empty, 1=nucleotide, 2=amino_acid"
    )`, It = `table defaultBedSchema
"BED12"
    (
    string chrom;      "The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)."
    uint   chromStart; "The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0."
    uint   chromEnd;   "The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99."
    string   name;   "Defines the name of the BED line."
    float   score;   "Feature score, doesn't care about the 0-1000 limit as in bed"
    char   strand;   "Defines the strand. Either '.' (=no strand) or '+' or '-'"
    uint thickStart; "The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position."
    uint thickEnd; "The ending position at which the feature is drawn thickly (for example the stop codon in gene displays)."
    string itemRgb; "An RGB value of the form R,G,B (e.g. 255,0,0). "
    uint blockCount; " The number of blocks (exons) in the BED line."
    uint[blockCount] blockSizes; " A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount."
    uint[blockCount] blockStarts; "A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount."
    )`, Gt = `table mafFrames
"codon frame assignment for MAF components"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start range in chromosome"
    uint   chromEnd;   "End range in chromosome"
    string src;        "Name of sequence source in MAF"
    ubyte frame;       "frame (0,1,2) for first base(+) or last bast(-)"
    char[1] strand;    "+ or -"
    string name;       "Name of gene used to define frame"
    int    prevFramePos;  "target position of the previous base (in transcription direction) that continues this frame, or -1 if none, or frame not contiguous"
    int    nextFramePos;  "target position of the next base (in transcription direction) that continues this frame, or -1 if none, or frame not contiguous"
    ubyte  isExonStart;  "does this start the CDS portion of an exon?"
    ubyte  isExonEnd;    "does this end the CDS portion of an exon?"
    )`, Ut = `table mafSummary
"Positions and scores for alignment blocks"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string src;        "Sequence name or database of alignment"
    float  score;      "Floating point score."
    char[1] leftStatus;  "Gap/break annotation for preceding block"
    char[1] rightStatus; "Gap/break annotation for following block"
    )`, jt = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  bigChain: Ft,
  bigGenePred: Rt,
  bigInteract: Tt,
  bigLink: zt,
  bigMaf: Bt,
  bigNarrowPeak: _t,
  bigPsl: Pt,
  defaultBedSchema: It,
  mafFrames: Gt,
  mafSummary: Ut
}, Symbol.toStringTag, { value: "Module" })), oe = Object.fromEntries(Object.entries(jt).map(([s, m]) => [
  s,
  Ge.parse(m.trim())
]));
function se(s) {
  const m = ["uint", "int", "float", "long"];
  return {
    ...s,
    fields: s.fields.map((t) => ({
      ...t,
      isArray: t.size && t.type !== "char",
      arrayIsNumeric: t.size && m.includes(t.type),
      isNumeric: !t.size && m.includes(t.type)
    }))
  };
}
const Mt = { ".": 0, "-": -1, "+": 1 };
function Ot(s) {
  var m;
  return s.length >= 12 && !Number.isNaN(parseInt(s[9], 10)) && ((m = s[10]) === null || m === void 0 ? void 0 : m.split(",").filter((t) => !!t).length) === parseInt(s[9], 10);
}
class Qt {
  constructor(m = {}) {
    if (m.autoSql)
      this.autoSql = se(Dt.parse(m.autoSql));
    else if (m.type) {
      if (!oe[m.type])
        throw new Error("Type not found");
      this.autoSql = se(oe[m.type]);
    } else
      this.autoSql = se(oe.defaultBedSchema), this.attemptDefaultBed = !0;
  }
  /*
   * parses a line of text as a BED line with the loaded autoSql schema
   *
   * @param line - a BED line as tab delimited text or array
   * @param opts - supply opts.uniqueId
   * @return a object representing a feature
   */
  parseLine(m, t = {}) {
    const { autoSql: d } = this, { uniqueId: w } = t, E = Array.isArray(m) ? m : m.split("	");
    let h = {};
    if (!this.attemptDefaultBed || this.attemptDefaultBed && Ot(E))
      for (let k = 0; k < d.fields.length; k++) {
        const y = d.fields[k];
        let a = E[k];
        const { isNumeric: f, isArray: p, arrayIsNumeric: v, name: A } = y;
        if (a == null)
          break;
        if (a !== ".") {
          if (f) {
            const q = Number(a);
            a = Number.isNaN(q) ? a : q;
          } else
            p && (a = a.split(","), a[a.length - 1] === "" && a.pop(), v && (a = a.map((q) => Number(q))));
          h[A] = a;
        }
      }
    else {
      const k = ["chrom", "chromStart", "chromEnd", "name"];
      h = Object.fromEntries(E.map((y, a) => [k[a] || "field" + a, y])), h.chromStart = +h.chromStart, h.chromEnd = +h.chromEnd, Number.isNaN(Number.parseFloat(h.field4)) || (h.score = +h.field4, delete h.field4), (h.field5 === "+" || h.field5 === "-") && (h.strand = h.field5, delete h.field5);
    }
    return w && (h.uniqueId = w), h.strand = Mt[h.strand] || 0, h.chrom = decodeURIComponent(h.chrom), h;
  }
}
function Lt(s) {
  const m = s.children(), t = s.get("thickStart"), d = s.get("thickEnd");
  if (!t && !d)
    return s;
  const w = m ? m.filter((y) => y.get("type") === "block").sort((y, a) => y.get("start") - a.get("start")) : [], E = [];
  w.forEach((y) => {
    const a = y.get("start"), f = y.get("end");
    if (t >= f) {
      const p = s.get("strand") > 0 ? "five" : "three";
      E.push({
        type: `${p}_prime_UTR`,
        start: a,
        end: f
      });
    } else if (t > a && t < f && d >= f) {
      const p = s.get("strand") > 0 ? "five" : "three";
      E.push({
        type: `${p}_prime_UTR`,
        start: a,
        end: t
      }, {
        type: "CDS",
        start: t,
        end: f
      });
    } else if (t <= a && d >= f)
      E.push({
        type: "CDS",
        start: a,
        end: f
      });
    else if (t > a && t < f && d < f) {
      const p = s.get("strand") > 0 ? "five" : "three", v = s.get("strand") > 0 ? "three" : "five";
      E.push({
        type: `${p}_prime_UTR`,
        start: a,
        end: t
      }, {
        type: "CDS",
        start: t,
        end: d
      }, {
        type: `${v}_prime_UTR`,
        start: d,
        end: f
      });
    } else if (t <= a && d > a && d < f) {
      const p = s.get("strand") > 0 ? "three" : "five";
      E.push({
        type: "CDS",
        start: a,
        end: d
      }, {
        type: `${p}_prime_UTR`,
        start: d,
        end: f
      });
    } else if (d <= a) {
      const p = s.get("strand") > 0 ? "three" : "five";
      E.push({
        type: `${p}_prime_UTR`,
        start: a,
        end: f
      });
    }
  });
  const h = Object.fromEntries(s.tags().map((y) => [y, s.get(y)]));
  return h.subfeatures = E, h.type = "mRNA", h.uniqueId = s.id(), delete h.chromStarts, delete h.chromStart, delete h.chromEnd, delete h.chrom, delete h.blockStarts, delete h.blockSizes, delete h.blockCount, delete h.thickStart, delete h.thickEnd, new Ie.SimpleFeature({
    data: h,
    id: s.id()
  });
}
function Vt(s, m) {
  return Object.fromEntries(m.split("	").map((t, d) => [s[d], t]));
}
function Wt(s, m, t, d, w, E, h, k) {
  const y = s.split("	"), a = y[m], f = +y[t], p = t === d ? 1 : 0, v = +y[d] + p, A = k ? Vt(k, s) : E.parseLine(s, { uniqueId: h }), { blockCount: q, blockSizes: X, blockStarts: I, chromStarts: O } = A;
  if (q) {
    const Y = O || I || [], L = X, H = f;
    A.subfeatures = [];
    for (let z = 0; z < q; z += 1) {
      const V = (Y[z] || 0) + H, U = V + (L[z] || 0);
      A.subfeatures.push({
        uniqueId: `${h}-${z}`,
        start: V,
        end: U,
        type: "block"
      });
    }
  }
  w && (A.score = +A[w]), delete A.chrom, delete A.chromStart, delete A.chromEnd;
  const G = new Ie.SimpleFeature({
    ...A,
    start: f,
    end: v,
    refName: a,
    uniqueId: h
  });
  return G.get("thickStart") ? Lt(G) : G;
}
export {
  Qt as B,
  Wt as f,
  Lt as u
};
