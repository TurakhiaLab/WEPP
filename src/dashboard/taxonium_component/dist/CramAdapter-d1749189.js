import { b3 as _, b4 as fr, aJ as cr, Q as Bt, b5 as Rt, b6 as lr, b7 as dr, b8 as hr, aV as ur, aZ as ht, aX as gr, aY as vr, aW as mr, t as wr } from "./index-8ee3baa4.js";
import { Q as Et } from "./index-05a6222b.js";
import { P as E } from "./binary_parser-421d5da5.js";
import { L as yr } from "./long-a46ae85c.js";
import { a as pr } from "./url-3d3a4b4f.js";
import "react";
import "react-dom";
function xt(n) {
  return _.Buffer.from(fr.inflate(n));
}
class nt extends Error {
}
class G extends Error {
}
class C extends nt {
}
class _r extends nt {
}
class Sr extends nt {
}
const x = 12, oe = 1 << x, Ee = 1 << 23;
class It {
  // int F, C;
  constructor() {
    this.F = void 0, this.C = void 0;
  }
}
class Ar {
  // final FC[] fc = new FC[256];
  // byte[] R;
  constructor() {
    this.fc = new Array(256);
    for (let e = 0; e < this.fc.length; e += 1)
      this.fc[e] = new It();
    this.R = null;
  }
}
let br = class {
  // int start; // Start of range.
  // int freq; // Symbol frequency.
  constructor() {
    this.start = void 0, this.freq = void 0;
  }
};
function Cr(n, e, t) {
  if (!(e <= 65536))
    throw new C("assertion failed: start <= 1<<16");
  if (!(t <= 65536 - e))
    throw new C("assertion failed: freq <= 1<<16");
  n.start = e, n.freq = t;
}
function Mt(n, e, t, r) {
  const a = (1 << r) - 1;
  return t * (n >> r) + (n & a) - e;
}
function Br(n, e, t) {
  return Mt(n, e.start, e.freq, t);
}
function Rr(n, e) {
  return n & (1 << e) - 1;
}
function Er(n, e, t, r, a) {
  const i = (1 << a) - 1;
  if (n = r * (n >> a) + (n & i) - t, n < Ee)
    do {
      const o = 255 & e.get();
      n = n << 8 | o;
    } while (n < Ee);
  return n;
}
function xr(n, e, t, r) {
  return Er(n, e, t.start, t.freq, r);
}
function Ir(n, e) {
  if (n < Ee)
    do
      n = n << 8 | 255 & e.get();
    while (n < Ee);
  return n;
}
const S = {
  FC: It,
  AriDecoder: Ar,
  Symbol: br,
  symbolInit: Cr,
  advanceStep: Mt,
  advanceSymbolStep: Br,
  get: Rr,
  advanceSymbol: xr,
  renormalize: Ir
};
function Ft(n) {
  if (!n)
    throw new C("assertion failed");
}
function Mr(n, e, t) {
  let r = 0, a = 0, i = n.get() & 255;
  do
    e.fc[i] == null && (e.fc[i] = new S.FC()), e.fc[i].F = n.get() & 255, e.fc[i].F >= 128 && (e.fc[i].F &= -129, e.fc[i].F = (e.fc[i].F & 127) << 8 | n.get() & 255), e.fc[i].C = a, S.symbolInit(t[i], e.fc[i].C, e.fc[i].F), e.R || (e.R = new Array(oe)), e.R.fill(i, a, a + e.fc[i].F), a += e.fc[i].F, r === 0 && i + 1 === (255 & n.getByteAt(n.position())) ? (i = n.get() & 255, r = n.get() & 255) : r !== 0 ? (r -= 1, i += 1) : i = n.get() & 255;
  while (i !== 0);
  Ft(a < oe);
}
function Fr(n, e, t) {
  let r = 0, a = 255 & n.get();
  do {
    let i = 0, o = 0, s = 255 & n.get();
    e[a] == null && (e[a] = new S.AriDecoder());
    do
      e[a].fc[s] == null && (e[a].fc[s] = new S.FC()), e[a].fc[s].F = 255 & n.get(), e[a].fc[s].F >= 128 && (e[a].fc[s].F &= -129, e[a].fc[s].F = (e[a].fc[s].F & 127) << 8 | 255 & n.get()), e[a].fc[s].C = o, e[a].fc[s].F === 0 && (e[a].fc[s].F = oe), t[a][s] == null && (t[a][s] = new S.RansDecSymbol()), S.symbolInit(t[a][s], e[a].fc[s].C, e[a].fc[s].F), e[a].R == null && (e[a].R = new Array(oe)), e[a].R.fill(s, o, o + e[a].fc[s].F), o += e[a].fc[s].F, Ft(o <= oe), i === 0 && s + 1 === (255 & n.getByteAt(n.position())) ? (s = 255 & n.get(), i = 255 & n.get()) : i !== 0 ? (i -= 1, s += 1) : s = 255 & n.get();
    while (s !== 0);
    r === 0 && a + 1 === (255 & n.getByteAt(n.position())) ? (a = 255 & n.get(), r = 255 & n.get()) : r !== 0 ? (r -= 1, a += 1) : a = 255 & n.get();
  } while (a !== 0);
}
function qr(n, e, t, r) {
  let a = n.getInt(), i = n.getInt(), o = n.getInt(), s = n.getInt();
  const c = r.remaining(), l = c & -4;
  for (let d = 0; d < l; d += 4) {
    const h = e.R[S.get(a, x)], u = e.R[S.get(i, x)], g = e.R[S.get(o, x)], m = e.R[S.get(s, x)];
    r.putAt(d, h), r.putAt(d + 1, u), r.putAt(d + 2, g), r.putAt(d + 3, m), a = S.advanceSymbolStep(a, t[255 & h], x), i = S.advanceSymbolStep(i, t[255 & u], x), o = S.advanceSymbolStep(o, t[255 & g], x), s = S.advanceSymbolStep(s, t[255 & m], x), a = S.renormalize(a, n), i = S.renormalize(i, n), o = S.renormalize(o, n), s = S.renormalize(s, n);
  }
  r.setPosition(l);
  let f;
  switch (c & 3) {
    case 0:
      break;
    case 1:
      f = e.R[S.get(a, x)], S.advanceSymbol(a, n, t[255 & f], x), r.put(f);
      break;
    case 2:
      f = e.R[S.get(a, x)], S.advanceSymbol(a, n, t[255 & f], x), r.put(f), f = e.R[S.get(i, x)], S.advanceSymbol(i, n, t[255 & f], x), r.put(f);
      break;
    case 3:
      f = e.R[S.get(a, x)], S.advanceSymbol(a, n, t[255 & f], x), r.put(f), f = e.R[S.get(i, x)], S.advanceSymbol(i, n, t[255 & f], x), r.put(f), f = e.R[S.get(o, x)], S.advanceSymbol(o, n, t[255 & f], x), r.put(f);
      break;
    default:
      throw new C("invalid output size encountered during rANS decoding");
  }
  r.setPosition(0);
}
function Pr(n, e, t, r) {
  const a = e.remaining();
  let i = n.getInt(), o = n.getInt(), s = n.getInt(), c = n.getInt();
  const l = a >> 2;
  let f = 0, d = l, h = 2 * l, u = 3 * l, g = 0, m = 0, v = 0, w = 0;
  for (; f < l; f += 1, d += 1, h += 1, u += 1) {
    const A = 255 & t[g].R[S.get(i, x)], M = 255 & t[m].R[S.get(o, x)], y = 255 & t[v].R[S.get(s, x)], p = 255 & t[w].R[S.get(c, x)];
    e.putAt(f, A), e.putAt(d, M), e.putAt(h, y), e.putAt(u, p), i = S.advanceSymbolStep(i, r[g][A], x), o = S.advanceSymbolStep(o, r[m][M], x), s = S.advanceSymbolStep(s, r[v][y], x), c = S.advanceSymbolStep(c, r[w][p], x), i = S.renormalize(i, n), o = S.renormalize(o, n), s = S.renormalize(s, n), c = S.renormalize(c, n), g = A, m = M, v = y, w = p;
  }
  for (; u < a; u += 1) {
    const A = 255 & t[w].R[S.get(c, x)];
    e.putAt(u, A), c = S.advanceSymbol(c, n, r[w][A], x), w = A;
  }
}
const Lr = 4;
function kr(n, e) {
  const t = new S.AriDecoder(), r = new Array(256);
  for (let a = 0; a < r.length; a += 1)
    r[a] = new S.Symbol();
  return Mr(n, t, r), qr(n, t, r, e), e;
}
function Ur(n, e) {
  const t = new Array(256);
  for (let a = 0; a < t.length; a += 1)
    t[a] = new S.AriDecoder();
  const r = new Array(256);
  for (let a = 0; a < r.length; a += 1) {
    r[a] = new Array(256);
    for (let i = 0; i < r[a].length; i += 1)
      r[a][i] = new S.Symbol();
  }
  return Fr(n, t, r), Pr(n, e, t, r), e;
}
class ut {
  constructor(e, t = 0) {
    this._buffer = e, this._position = t, this.length = e.length;
  }
  get() {
    const e = this._buffer[this._position];
    return this._position += 1, e;
  }
  getByte() {
    return this.get();
  }
  getByteAt(e) {
    return this._buffer[e];
  }
  position() {
    return this._position;
  }
  put(e) {
    return this._buffer[this._position] = e, this._position += 1, e;
  }
  putAt(e, t) {
    return this._buffer[e] = t, t;
  }
  setPosition(e) {
    return this._position = e, e;
  }
  getInt() {
    const e = this._buffer.readInt32LE(this._position);
    return this._position += 4, e;
  }
  remaining() {
    return this._buffer.length - this._position;
  }
}
function Nr(n, e, t = 0) {
  if (n.length === 0)
    return e.fill(0), e;
  const r = new ut(n, t), a = r.get();
  if (a !== 0 && a !== 1)
    throw new C(`Invalid rANS order ${a}`);
  if (r.getInt() !== r.remaining() - Lr)
    throw new C("Incorrect input length.");
  const o = r.getInt(), s = new ut(e || _.Buffer.allocUnsafe(o));
  if (s.length < o)
    throw new C(`Output buffer too small to fit ${o} bytes.`);
  switch (a) {
    case 0:
      return kr(r, s);
    case 1:
      return Ur(r, s);
    default:
      throw new C(`Invalid rANS order: ${a}`);
  }
}
const we = new E().itf8(), be = {
  parser: new E().string("magic", { length: 4 }).uint8("majorVersion").uint8("minorVersion").string("fileId", { length: 20, stripNull: !0 }),
  maxLength: 26
}, Dr = {
  parser: new E().uint8("compressionMethod", {
    formatter: (
      /* istanbul ignore next */
      /* istanbul ignore next */
      (n) => {
        const e = [
          "raw",
          "gzip",
          "bzip2",
          "lzma",
          "rans",
          "rans4x16",
          "arith",
          "fqzcomp",
          "tok3"
        ][n];
        if (!e)
          throw new Error(`compression method number ${n} not implemented`);
        return e;
      }
    )
  }).uint8("contentType", {
    formatter: (
      /* istanbul ignore next */
      /* istanbul ignore next */
      (n) => {
        const e = [
          "FILE_HEADER",
          "COMPRESSION_HEADER",
          "MAPPED_SLICE_HEADER",
          "UNMAPPED_SLICE_HEADER",
          "EXTERNAL_DATA",
          "CORE_DATA"
        ][n];
        if (!e)
          throw new Error(`invalid block content type id ${n}`);
        return e;
      }
    )
  }).itf8("contentId").itf8("compressedSize").itf8("uncompressedSize"),
  maxLength: 17
}, $r = {
  parser: new E().uint32("crc32"),
  maxLength: 4
}, zr = new E().itf8("size").buffer("ents", {
  length: "size",
  formatter: (
    /* istanbul ignore next */
    /* istanbul ignore next */
    (n) => {
      function e(i, o) {
        const s = n.toString("utf8", i, o), c = [];
        for (let l = 0; l < s.length; l += 3)
          c.push(s.substr(l, 3));
        return c;
      }
      var t = [], r = 0, a;
      for (a = 0; a < n.length; a += 1)
        n[a] || (t.push(e(r, a)), r = a + 1);
      return a > r && t.push(e(r, a)), t;
    }
  )
}), j = new E().uint8(null, {
  formatter: (
    /* istanbul ignore next */
    /* istanbul ignore next */
    (n) => !!n
  )
}), Wr = new E().itf8("mapSize").itf8("mapCount").array("ents", {
  length: "mapCount",
  type: new E().string("key", {
    length: 2,
    stripNull: !1
    // formatter: val => cramPreservationMapKeys[val] || 0,
  }).choice("value", {
    tag: "key",
    choices: {
      MI: j,
      UI: j,
      PI: j,
      RN: j,
      AP: j,
      RR: j,
      SM: new E().array(null, { type: "uint8", length: 5 }),
      TD: new E().nest(null, {
        type: zr,
        formatter: (
          /* istanbul ignore next */
          /* istanbul ignore next */
          (n) => n.ents
        )
      })
    }
  })
});
function Ue(n) {
  const e = {};
  for (let t = 0; t < n.ents.length; t += 1) {
    const { key: r, value: a } = n.ents[t];
    e[r] && console.warn(`duplicate key ${r} in map`), e[r] = a;
  }
  return e;
}
const Tr = {
  cramFileDefinition: be,
  cramBlockHeader: Dr,
  cramBlockCrc32: $r
};
function ie(n) {
  return typeof n.refSeqId == "number";
}
const gt = {
  // assemble a section parser for the unmapped slice header, with slight
  // variations depending on the major version of the cram file
  cramUnmappedSliceHeader(n) {
    let e = 0, t = new E().itf8("numRecords");
    return e += 5, n >= 3 ? (t = t.ltf8("recordCounter"), e += 9) : n === 2 && (t = t.itf8("recordCounter"), e += 5), t = t.itf8("numBlocks").itf8("numContentIds").array("contentIds", {
      type: we,
      length: "numContentIds"
    }), e += 5 * 2, n >= 2 && (t = t.array("md5", { type: "uint8", length: 16 }), e += 16), { parser: t, maxLength: (a) => e + a * 5 };
  },
  // assembles a section parser for the unmapped slice header, with slight
  // variations depending on the major version of the cram file
  cramMappedSliceHeader(n) {
    let e = new E().itf8("refSeqId").itf8("refSeqStart").itf8("refSeqSpan").itf8("numRecords"), t = 5 * 4;
    return n >= 3 ? (e = e.ltf8("recordCounter"), t += 9) : n === 2 && (e = e.itf8("recordCounter"), t += 5), e = e.itf8("numBlocks").itf8("numContentIds").array("contentIds", {
      type: we,
      length: "numContentIds"
    }).itf8("refBaseBlockId"), t += 5 * 3, n >= 2 && (e = e.array("md5", { type: "uint8", length: 16 }), t += 16), { parser: e, maxLength: (a) => t + a * 5 };
  },
  cramEncoding(n) {
    return { parser: new E().namely("cramEncoding").itf8("codecId").itf8("parametersBytes").choice("parameters", {
      tag: "codecId",
      choices: {
        0: new E(),
        1: new E().itf8("blockContentId"),
        2: new E().itf8("offset").itf8("M"),
        // HUFFMAN_INT
        3: E.start().itf8("numCodes").array("symbols", { length: "numCodes", type: we }).itf8("numLengths").array("bitLengths", { length: "numLengths", type: we }),
        4: E.start().nest("lengthsEncoding", { type: "cramEncoding" }).nest("valuesEncoding", { type: "cramEncoding" }),
        // BYTE_ARRAY_STOP is a little different for CRAM v1
        5: new E().uint8("stopByte").itf8("blockContentId"),
        6: new E().itf8("offset").itf8("length"),
        7: new E().itf8("offset").itf8("K"),
        8: new E().itf8("offset").itf8("log2m"),
        9: new E().itf8("offset")
        // GAMMA
      }
    }) };
  },
  cramDataSeriesEncodingMap(n) {
    return new E().itf8("mapSize").itf8("mapCount").array("ents", {
      length: "mapCount",
      type: new E().string("key", { length: 2, stripNull: !1 }).nest("value", { type: this.cramEncoding(n).parser })
    });
  },
  cramTagEncodingMap(n) {
    return new E().itf8("mapSize").itf8("mapCount").array("ents", {
      length: "mapCount",
      type: new E().itf8("key", {
        formatter: (
          /* istanbul ignore next */
          /* istanbul ignore next */
          (e) => (
            /* istanbul ignore next */
            String.fromCharCode(e >> 16 & 255) + String.fromCharCode(e >> 8 & 255) + String.fromCharCode(e & 255)
          )
        )
      }).nest("value", { type: this.cramEncoding(n).parser })
    });
  },
  cramCompressionHeader(n) {
    let e = new E();
    return e = e.nest("preservation", {
      type: Wr,
      formatter: Ue
    }).nest("dataSeriesEncoding", {
      type: this.cramDataSeriesEncodingMap(n),
      formatter: Ue
    }).nest("tagEncoding", {
      type: this.cramTagEncodingMap(n),
      formatter: Ue
    }), { parser: e };
  },
  cramContainerHeader1(n) {
    let e = new E().int32("length").itf8("refSeqId").itf8("refSeqStart").itf8("alignmentSpan").itf8("numRecords"), t = 4 + 5 * 4;
    return n >= 3 ? (e = e.ltf8("recordCounter"), t += 9) : n === 2 && (e = e.itf8("recordCounter"), t += 5), n > 1 && (e = e.ltf8("numBases"), t += 9), e = e.itf8("numBlocks").itf8("numLandmarks"), t += 5 + 5, { parser: e, maxLength: t };
  },
  cramContainerHeader2(n) {
    let e = new E().itf8("numLandmarks").array("landmarks", {
      type: new E().itf8(),
      length: "numLandmarks"
    }), t = 0;
    return n >= 3 && (e = e.uint32("crc32"), t = 4), {
      parser: e,
      maxLength: (r) => 5 + r * 5 + t
    };
  }
};
function Hr(n) {
  const e = Object.assign({}, Tr);
  return Object.keys(gt).forEach((t) => {
    e[t] = gt[t](n);
  }), e;
}
var ve = class {
  constructor(e, t = 0, r = 0) {
    r != 0 ? (this.buf = _.Buffer.allocUnsafe(r), this.length = r) : (this.buf = e, this.length = e.length), this.pos = t;
  }
  // ----------
  // Reading
  EOF() {
    return this.pos >= this.length;
  }
  ReadData(e) {
    var t = this.buf.slice(this.pos, this.pos + e);
    return this.pos += e, t;
  }
  ReadByte() {
    const e = this.buf[this.pos];
    return this.pos++, e;
  }
  ReadChar() {
    const e = this.buf[this.pos];
    return this.pos++, String.fromCharCode(e);
  }
  ReadUint16() {
    var e = this.ReadByte();
    return e |= this.ReadByte() << 8, e;
  }
  ReadUint32() {
    const e = this.buf.readInt32LE(this.pos);
    return this.pos += 4, e;
  }
  // nul terminated string
  ReadString() {
    var e = "";
    do {
      var t = this.buf[this.pos++];
      t && (e += String.fromCharCode(t));
    } while (t);
    return e;
  }
  //    ReadUint7() {
  //	// Variable sized unsigned integers
  //	var i = 0;
  //	var s = 0;
  //	do {
  //	    var c = this.ReadByte();
  //	    i = i | ((c & 0x7f)<<s);
  //	    s += 7;
  //	} while ((c & 0x80))
  //
  //	return i;
  //    }
  ReadUint7() {
    var e = 0;
    do {
      var t = this.ReadByte();
      e = e << 7 | t & 127;
    } while (t & 128);
    return e;
  }
  ReadITF8() {
    var e = this.buf[this.pos];
    return this.pos++, e >= 240 ? (e = (e & 15) << 28, e += (this.buf[this.pos + 0] << 20) + (this.buf[this.pos + 1] << 12) + (this.buf[this.pos + 2] << 4) + (this.buf[this.pos + 3] >> 4), this.pos += 4) : e >= 224 ? (e = (e & 15) << 24, e += (this.buf[this.pos + 0] << 16) + (this.buf[this.pos + 1] << 8) + (this.buf[this.pos + 2] << 0), this.pos += 3) : e >= 192 ? (e = (e & 31) << 16, e += (this.buf[this.pos + 0] << 8) + (this.buf[this.pos + 1] << 0), this.pos += 2) : e >= 128 && (e = (e & 63) << 8, e += this.buf[this.pos], this.pos++), e;
  }
  // ----------
  // Writing
  WriteByte(e) {
    this.buf[this.pos++] = e;
  }
  WriteChar(e) {
    this.buf[this.pos++] = e.charCodeAt(0);
  }
  WriteString(e) {
    for (var t = 0; t < e.length; t++)
      this.buf[this.pos++] = e.charCodeAt(t);
    this.buf[this.pos++] = 0;
  }
  WriteData(e, t) {
    for (var r = 0; r < t; r++)
      this.buf[this.pos++] = e[r];
  }
  WriteStream(e) {
    this.WriteData(e.buf, e.pos);
  }
  WriteUint16(e) {
    this.WriteByte(e & 255), this.WriteByte(e >> 8 & 255);
  }
  WriteUint32(e) {
    this.buf.writeInt32LE(e, this.pos), this.pos += 4;
  }
  //    WriteUint7(i) {
  //	do {
  //	    this.WriteByte((i & 0x7f) | ((i > 0x80) << 7));
  //	    i >>= 7;
  //	} while (i > 0);
  //    }
  WriteUint7(e) {
    var t = 0, r = e;
    do
      t += 7, r >>= 7;
    while (r > 0);
    do
      t -= 7, this.WriteByte((e >> t & 127) + ((t > 0) << 7));
    while (t > 0);
  }
  WriteITF8(e) {
    e < 0 && (e = 1 + e), e <= 127 ? this.buf[this.pos++] = e : e <= 16383 ? (this.buf[this.pos++] = 128 | Math.floor(e / 256), this.buf[this.pos++] = e & 255) : e < 131071 ? (this.buf[this.pos++] = 192 | Math.floor(e / 65536), this.buf[this.pos++] = Math.floor(e / 256) & 255, this.buf[this.pos++] = e & 255) : e < 268435455 ? (this.buf[this.pos++] = 224 | Math.floor(e / 16777216), this.buf[this.pos++] = Math.floor(e / 65536) & 255, this.buf[this.pos++] = Math.floor(e / 256) & 255, this.buf[this.pos++] = e & 255) : (this.buf[this.pos++] = 240 | Math.floor(e / 268435456), this.buf[this.pos++] = Math.floor(e / 1048576) & 255, this.buf[this.pos++] = Math.floor(e / 4096) & 255, this.buf[this.pos++] = Math.floor(e / 4) & 255, this.buf[this.pos++] = e & 15);
  }
  // ----------
  // Writing from end of buffer going backwards.
  // Needed by rANS codec.
  WriteByteNeg(e) {
    this.buf[--this.pos] = e;
  }
};
const fe = ve;
function $e(n) {
  return n & 4095;
}
function Or(n, e) {
  for (var t = 0; e >= n[t + 1]; )
    t++;
  return t;
}
function qt(n) {
  for (var e = new Array(4096), t = 0, r = 0; r < 4096; r++) {
    for (; r >= n[t + 1]; )
      t++;
    e[r] = t;
  }
  return e;
}
function ze(n, e, t) {
  return t * (n >> 12) + (n & 4095) - e;
}
function We(n, e) {
  for (; e < 1 << 23; )
    e = (e << 8) + n.ReadByte();
  return e;
}
function Pt() {
  return 1 << 23;
}
function Lt(n, e) {
  e.WriteByteNeg(n >> 24 & 255), e.WriteByteNeg(n >> 16 & 255), e.WriteByteNeg(n >> 8 & 255), e.WriteByteNeg(n >> 0 & 255);
}
function Gr(n, e, t, r) {
  for (var a = (8388608 >> r << 8) * t; n >= a; )
    e.WriteByteNeg(n & 255), n >>= 8;
  return n;
}
function Ce(n, e, t, r, a) {
  return n = Gr(n, e, r, a), n = (Math.floor(n / r) << a) + n % r + t, n;
}
function Qr(n) {
  var e = new fe(n), t = e.ReadByte();
  e.ReadUint32();
  var r = e.ReadUint32();
  return t == 0 ? Kr(e, r) : Jr(e, r);
}
function Vr(n, e) {
  return e == 0 ? Yr(n) : rn(n);
}
function kt(n, e, t) {
  for (var r = 0; r < 256; r++)
    e[r] = 0;
  var a = n.ReadByte(), i = a, o = 0;
  do {
    var s = n.ReadITF8();
    e[a] = s, o > 0 ? (o--, a++) : (a = n.ReadByte(), a == i + 1 && (o = n.ReadByte())), i = a;
  } while (a != 0);
  t[0] = 0;
  for (var r = 0; r <= 255; r++)
    t[r + 1] = t[r] + e[r];
}
function Kr(n, e) {
  var t = new Array(256), r = new Array(256);
  kt(n, t, r);
  for (var a = qt(r), i = new Array(4), o = 0; o < 4; o++)
    i[o] = n.ReadUint32();
  for (var s = new _.Buffer.allocUnsafe(e), o = 0; o < e; o++) {
    var c = o % 4, l = $e(i[c]), f = a[l];
    s[o] = f, i[c] = ze(i[c], r[f], t[f]), i[c] = We(n, i[c]);
  }
  return s;
}
function Xr(n, e) {
  for (var t = 0; t < 256; t++)
    e[t] = 0;
  for (var t = 0; t < n.length; t++)
    e[n[t]]++;
}
function Ut(n) {
  for (var e = 0, t = 0; t < 256; t++)
    e += n[t];
  const r = 4096;
  var a = r / e;
  do {
    var i = 0, o = 0, s = 0;
    e = 0;
    for (var t = 0; t < 256; t++)
      n[t] != 0 && (i < n[t] && (i = n[t], o = t), n[t] = Math.floor(n[t] * a), n[t] == 0 && (n[t] = 1), e += n[t]);
    e < r ? n[o] += r - e : e - r < n[o] / 2 && n[o] > 2 ? n[o] -= e - r : e != r && (a = a * 0.99, s = 1);
  } while (s);
}
function Nt(n, e) {
  for (var t = 0, r = 0; r < 256; r++)
    if (e[r]) {
      if (t > 0)
        t--;
      else if (n.WriteByte(r), r > 0 && e[r - 1] > 0) {
        for (t = r + 1; t < 256 && e[t]; t++)
          ;
        t -= r + 1, n.WriteByte(t);
      }
      n.WriteITF8(e[r]);
    }
  n.WriteByte(0);
}
function Yr(n) {
  const e = n.length;
  var t = new fe("", 0, 257 * 3 + 9);
  t.WriteByte(0), t.WriteUint32(0), t.WriteUint32(0);
  var r = new Array(256);
  Xr(n, r), Ut(r), Nt(t, r);
  var a = new Array(256);
  a[0] = 0;
  for (var i = 1; i < 256; i++)
    a[i] = a[i - 1] + r[i - 1];
  for (var o = new Array(4), i = 0; i < 4; i++)
    o[i] = Pt();
  for (var s = Math.floor(e * 1.05 + 100), c = new fe("", s, s), i = e - 1; i >= 0; i--)
    o[i % 4] = Ce(o[i % 4], c, a[n[i]], r[n[i]], 12);
  for (var i = 3; i >= 0; i--)
    Lt(o[i], c);
  var l = t.pos;
  return t.buf.writeInt32LE(l - 9 + (c.length - c.pos), 1), t.buf.writeInt32LE(e, 5), _.Buffer.concat(
    [
      t.buf.slice(0, t.pos),
      c.buf.slice(c.pos, c.length)
    ],
    t.pos + c.length - c.pos
  );
}
function Zr(n, e, t) {
  for (var r = 0; r < 256; r++) {
    e[r] = new Array(256), t[r] = new Array(256);
    for (var a = 0; a < 256; a++)
      e[r][a] = 0;
  }
  var i = n.ReadByte(), o = i, s = 0;
  do
    kt(n, e[i], t[i]), s > 0 ? (s--, i++) : (i = n.ReadByte(), i == o + 1 && (s = n.ReadByte())), o = i;
  while (i != 0);
}
function Jr(n, e) {
  var t = new Array(256), r = new Array(256);
  Zr(n, t, r);
  for (var a = new Array(256), i = 0; i < 256; i++)
    a[i] = qt(r[i]);
  for (var o = new Array(4), s = new Array(4), c = 0; c < 4; c++)
    o[c] = n.ReadUint32(), s[c] = 0;
  for (var l = new _.Buffer.allocUnsafe(e), f = Math.floor(e / 4), i = 0; i < f; i++)
    for (var c = 0; c < 4; c++) {
      var d = $e(o[c]), h = a[s[c]][d];
      l[i + c * f] = h, o[c] = ze(o[c], r[s[c]][h], t[s[c]][h]), o[c] = We(n, o[c]), s[c] = h;
    }
  for (i = 4 * i; i < e; ) {
    var d = $e(o[3]), h = Or(r[s[3]], d);
    l[i++] = h, o[3] = ze(o[3], r[s[3]][h], t[s[3]][h]), o[3] = We(n, o[3]), s[3] = h;
  }
  return l;
}
function jr(n, e, t) {
  for (var r = 0; r < 256; r++) {
    t[r] = 0;
    for (var a = 0; a < 256; a++)
      e[r][a] = 0;
  }
  for (var i = 0, r = 0; r < n.length; r++)
    t[n[r]]++, e[i][n[r]]++, i = n[r];
  e[0][n[1 * (n.length >> 2)]]++, e[0][n[2 * (n.length >> 2)]]++, e[0][n[3 * (n.length >> 2)]]++, t[0] += 3;
}
function en(n, e) {
  for (var t = 0; t < 256; t++)
    e[t] && Ut(n[t]);
}
function tn(n, e, t) {
  for (var r = 0, a = 0; a < 256; a++)
    if (t[a]) {
      if (r > 0)
        r--;
      else if (n.WriteByte(a), a > 0 && t[a - 1] > 0) {
        for (r = a + 1; r < 256 && t[r]; r++)
          ;
        r -= a + 1, n.WriteByte(r);
      }
      Nt(n, e[a]);
    }
  n.WriteByte(0);
}
function rn(n) {
  const e = n.length;
  var t = new fe("", 0, 257 * 257 * 3 + 9);
  t.WriteByte(1), t.WriteUint32(0), t.WriteUint32(0);
  for (var r = new Array(256), a = new Array(256), i = new Array(256), o = 0; o < 256; o++)
    a[o] = new Array(256), i[o] = new Array(256);
  jr(n, a, r), en(a, r), tn(t, a, r);
  for (var o = 0; o < 256; o++)
    if (r[o]) {
      i[o][0] = 0;
      for (var s = 1; s < 256; s++)
        i[o][s] = i[o][s - 1] + a[o][s - 1];
    }
  for (var c = new Array(4), l = new Array(4), s = 0; s < 4; s++)
    c[s] = Pt(), l[s] = 0;
  for (var f = new fe("", e, e), d = Math.floor(e / 4), h = new Array(4), u = new Array(4), s = 0; s < 4; s++)
    h[s] = (s + 1) * d - 2, u[s] = n[h[s] + 1];
  u[3] = n[e - 1];
  for (var o = e - 2; o > 4 * d - 2; o--)
    c[3] = Ce(c[3], f, i[n[o]][u[3]], a[n[o]][u[3]], 12), u[3] = n[o];
  for (; h[0] >= 0; )
    for (var s = 3; s >= 0; s--) {
      var g = n[h[s]];
      c[s] = Ce(c[s], f, i[g][u[s]], a[g][u[s]], 12), u[s] = g, h[s]--;
    }
  for (var s = 3; s >= 0; s--)
    c[s] = Ce(c[s], f, i[0][u[s]], a[0][u[s]], 12);
  for (var o = 3; o >= 0; o--)
    Lt(c[o], f);
  var m = t.pos;
  return t.buf.writeInt32LE(m - 9 + (f.length - f.pos), 1), t.buf.writeInt32LE(e, 5), _.Buffer.concat(
    [
      t.buf.slice(0, t.pos),
      f.buf.slice(f.pos, f.length)
    ],
    t.pos + f.length - f.pos
  );
}
var nn = { decode: Qr, encode: Vr };
const U = ve;
function Te(n, e) {
  return n & (1 << e) - 1;
}
function an(n, e) {
  for (var t = 0; e >= n[t + 1]; )
    t++;
  return t;
}
function Dt(n, e) {
  for (var t = 1 << e, r = new Array(t), a = 0, i = 0; i < t; i++) {
    for (; i >= n[a + 1]; )
      a++;
    r[i] = a;
  }
  return r;
}
function He(n, e, t, r) {
  return t * (n >> r) + (n & (1 << r) - 1) - e;
}
function Oe(n, e) {
  return e < 32768 && (e = (e << 16) + n.ReadUint16()), e;
}
function $t() {
  return 32768;
}
function zt(n, e) {
  e.WriteByteNeg(n >> 24 & 255), e.WriteByteNeg(n >> 16 & 255), e.WriteByteNeg(n >> 8 & 255), e.WriteByteNeg(n >> 0 & 255);
}
function on(n, e, t, r) {
  for (var a = (1 << 31 - r) * t; n >= a; )
    e.WriteByteNeg(n >> 8 & 255), e.WriteByteNeg(n & 255), n >>= 16;
  return n;
}
function Be(n, e, t, r, a) {
  return n = on(n, e, r, a), n = (Math.floor(n / r) << a) + n % r + t, n;
}
function sn(n) {
  for (var e = new Array(256), t = 0; t < 256; t++)
    e[t] = 0;
  for (var r = -1, t = 0; t < n.length; t++)
    e[n[t]] += n[t] == r ? 1 : -1, r = n[t];
  for (var a = 0, t = 0; t < 256; t++)
    e[t] > 0 && a++;
  a || (a = 1, e[0] = 1);
  var f = new U("", 0, a + 1 + n.length);
  f.WriteByte(a);
  for (var t = 0; t < 256; t++)
    e[t] > 0 && f.WriteByte(t);
  for (var i = new _.Buffer.allocUnsafe(n.length), o = 0, t = 0; t < n.length; t++)
    if (i[o++] = n[t], e[n[t]] > 0) {
      r = n[t];
      for (var s = 0; t + s + 1 < n.length && n[t + s + 1] == r; )
        s++;
      f.WriteUint7(s), t += s;
    }
  var c = ot(f.buf.slice(0, f.pos)), l = new U("", 0, 16);
  l.WriteUint7(f.pos * 2), l.WriteUint7(o), l.WriteUint7(c.length);
  var f = _.Buffer.concat([l.buf.slice(0, l.pos), c]);
  return [f, i.slice(0, o)];
}
function fn(n) {
  var e = n.ReadUint7(), t = n.ReadUint7();
  if (e & 1)
    var a = n.ReadData((e - 1) / 2);
  else {
    var r = n.ReadUint7(), a = n.ReadData(r);
    a = at(new U(a), e / 2);
  }
  var a = new U(a), i = new Array(256), o = a.ReadByte();
  o == 0 && (o = 256);
  for (var s = 0; s < o; s++)
    i[a.ReadByte()] = 1;
  return [i, a, t];
}
function cn(n, e, t, r) {
  new U(n);
  for (var a = new _.Buffer.allocUnsafe(r), i = 0, o = 0; i < r; o++) {
    var s = n[o];
    if (e[s])
      for (var c = t.ReadUint7(), l = 0; l <= c; l++)
        a[i++] = s;
    else
      a[i++] = s;
  }
  return a;
}
function ln(n) {
  for (var e = new Array(256), t = 0; t < 256; t++)
    e[t] = 0;
  for (var t = 0; t < n.length; t++)
    e[n[t]]++;
  for (var r = new Array(256), a = 0, t = 0; t < 256; t++)
    e[t] > 0 && (r[t] = a++);
  if (!(a > 16)) {
    if (a <= 1)
      var i = new _.Buffer.allocUnsafe(0);
    else if (a <= 2) {
      var i = new _.Buffer.allocUnsafe(Math.ceil(n.length / 8)), s = -1;
      for (t = 0; t < n.length; t++)
        t % 8 == 0 && (i[++s] = 0), i[s] += r[n[t]] << t % 8;
    } else if (a <= 4) {
      var i = new _.Buffer.allocUnsafe(Math.ceil(n.length / 4)), s = -1;
      for (t = 0; t < n.length; t++)
        t % 4 == 0 && (i[++s] = 0), i[s] += r[n[t]] << t % 4 * 2;
    } else {
      var i = new _.Buffer.allocUnsafe(Math.ceil(n.length / 2)), s = -1;
      for (t = 0; t < n.length; t++)
        t % 2 == 0 && (i[++s] = 0), i[s] += r[n[t]] << t % 2 * 4;
    }
    var o = new U("", 0, a + 5);
    o.WriteByte(a);
    for (var s = 0, t = 0; t < 256; t++)
      e[t] > 0 && (e[t] = s++, o.WriteByte(t));
    return o.WriteUint7(i.length), [o.buf.slice(0, o.pos), i];
  }
}
function dn(n) {
  for (var e = n.ReadByte(), t = new Array(e), r = 0; r < e; r++)
    t[r] = n.ReadByte();
  var a = n.ReadUint7();
  return [t, e, a];
}
function hn(n, e, t, r) {
  var a = new _.Buffer.allocUnsafe(r), i = 0;
  if (t <= 1)
    for (var o = 0; o < r; o++)
      a[o] = e[0];
  else if (t <= 2)
    for (o = 0; o < r; o++) {
      if (o % 8 == 0)
        var s = n[i++];
      a[o] = e[s & 1], s >>= 1;
    }
  else if (t <= 4)
    for (o = 0; o < r; o++) {
      if (o % 4 == 0)
        var s = n[i++];
      a[o] = e[s & 3], s >>= 2;
    }
  else if (t <= 16)
    for (o = 0; o < r; o++) {
      if (o % 2 == 0)
        var s = n[i++];
      a[o] = e[s & 15], s >>= 4;
    }
  return a;
}
function un(n, e, t) {
  t == 0 && (t = 4);
  for (var r = new Array(t), a = new Array(t), i = 0; i < t; i++)
    a[i] = Math.floor(e.length / t) + (e.length % t > i), r[i] = new Array(a[i]);
  for (var o = 0, s = 0; s < e.length; s += t, o++)
    for (var c = 0; c < t; c++)
      o < r[c].length && (r[c][o] = e[s + c]);
  for (var l = new Array(t), f = 0, i = 0; i < t; i++) {
    var d = Ge(r[i], 0), h = Ge(r[i], 1);
    l[i] = h.length < d.length ? h : d, f += l[i].length;
  }
  var u = new U("", 0, f + 5 * t + 1);
  u.WriteByte(t);
  for (var i = 0; i < t; i++)
    u.WriteUint7(l[i].length);
  for (var i = 0; i < t; i++)
    u.WriteData(l[i], l[i].length);
  return u.buf.slice(0, u.buf.pos);
}
function gn(n, e) {
  for (var t = n.ReadByte(), r = new Array(t), a = new Array(t), i = 0; i < t; i++)
    r[i] = n.ReadUint7();
  for (var o = new Array(t), i = 0; i < t; i++)
    a[i] = Math.floor(e / t) + (e % t > i), o[i] = Wt(n, a[i]);
  for (var s = new _.Buffer.allocUnsafe(e), i = 0; i < t; i++)
    for (var c = 0; c < a[i]; c++)
      s[c * t + i] = o[i][c];
  return s;
}
function vn(n) {
  var e = new U(n);
  return Wt(e, 0);
}
function Wt(n, e) {
  var t = n.ReadByte(), r = t & 1, a = t & 8, i = t & 16, o = t & 32, s = t & 64, c = t & 128;
  if (i || (e = n.ReadUint7()), a)
    return gn(n, e);
  if (c)
    var l = e, [f, d, e] = dn(n);
  if (s)
    var h = e, [u, g, e] = fn(n);
  if (o)
    var m = n.ReadData(e);
  else if (r == 0)
    var m = at(n, e);
  else
    var m = _n(n, e);
  return s && (m = cn(m, u, g, h)), c && (m = hn(m, f, d, l)), m;
}
function Ge(n, e) {
  var t = new U("", 0, 10);
  t.WriteByte(e);
  var r = e & 1, a = e & 8, i = e & 16, o = e & 32, s = e & 64, c = e & 128, l = e >> 8;
  if (i || t.WriteUint7(n.length), a)
    return _.Buffer.concat([t.buf.slice(0, t.pos), un(t, n, l)]);
  var f = new _.Buffer.alloc(0);
  c && ([f, n] = ln(n));
  var d = new _.Buffer.alloc(0);
  if (s && ([d, n] = sn(n)), n.length < 4 && r == 1 && (r = 0, t.buf[0] &= -2), o)
    var h = n;
  else if (r == 0)
    var h = ot(n);
  else
    var h = Bn(n);
  return _.Buffer.concat([t.buf.slice(0, t.pos), f, d, h]);
}
function Tt(n) {
  for (var e = new Array(256), t = 0; t < 256; t++)
    e[t] = 0;
  var r = 0, a = n.ReadByte(), i = a;
  do
    e[a] = 1, r > 0 ? (r--, a++) : (a = n.ReadByte(), a == i + 1 && (r = n.ReadByte())), i = a;
  while (a != 0);
  return e;
}
function mn(n, e, t) {
  for (var r = 0; r < 256; r++)
    e[r] = 0;
  for (var a = Tt(n), r = 0; r < 256; r++)
    a[r] > 0 && (e[r] = n.ReadUint7());
  it(e, 12), t[0] = 0;
  for (var r = 0; r <= 255; r++)
    t[r + 1] = t[r] + e[r];
}
function at(n, e) {
  var t = new Array(256), r = new Array(256);
  mn(n, t, r);
  for (var a = Dt(r, 12), i = new Array(4), o = 0; o < 4; o++)
    i[o] = n.ReadUint32();
  for (var s = new _.Buffer.allocUnsafe(e), o = 0; o < e; o++) {
    var c = o % 4, l = Te(i[c], 12), f = a[l];
    s[o] = f, i[c] = He(i[c], r[f], t[f], 12), i[c] = Oe(n, i[c]);
  }
  return s;
}
function wn(n, e) {
  for (var t = 0; t < 256; t++)
    e[t] = 0;
  for (var t = 0; t < n.length; t++)
    e[n[t]]++;
}
function Qe(n, e) {
  for (var t = 0, r = 0; r < 256; r++)
    t += n[r];
  const a = 1 << e;
  var i = a / t;
  do {
    var o = 0, s = 0, c = 0;
    t = 0;
    for (var r = 0; r < 256; r++)
      n[r] != 0 && (o < n[r] && (o = n[r], s = r), n[r] = Math.floor(n[r] * i), n[r] == 0 && (n[r] = 1), t += n[r]);
    t < a ? n[s] += a - t : t - a < n[s] / 2 && n[s] > 2 ? n[s] -= t - a : t != a && (i = a / t, c = 1);
  } while (c);
}
function it(n, e) {
  for (var t = 0, r = 0; r < 256; r++)
    t += n[r];
  if (!(t == 0 || t == 1 << e)) {
    for (var a = 0; t < 1 << e; )
      t *= 2, a++;
    for (var r = 0; r < 256; r++)
      n[r] <<= a;
  }
}
function Ht(n, e) {
  for (var t = 0, r = 0; r < 256; r++)
    if (e[r]) {
      if (t > 0)
        t--;
      else if (n.WriteByte(r), r > 0 && e[r - 1] > 0) {
        for (t = r + 1; t < 256 && e[t]; t++)
          ;
        t -= r + 1, n.WriteByte(t);
      }
    }
  n.WriteByte(0);
}
function yn(n, e) {
  Ht(n, e);
  for (var t = 0; t < 256; t++)
    e[t] && n.WriteUint7(e[t]);
}
function ot(n) {
  const e = n.length;
  var t = new U("", 0, 257 * 3 + 9), r = new Array(256);
  wn(n, r);
  var a = Math.ceil(Math.log2(e));
  a > 12 && (a = 12), Qe(r, a), yn(t, r), Qe(r, 12);
  var i = new Array(256);
  i[0] = 0;
  for (var o = 1; o < 256; o++)
    i[o] = i[o - 1] + r[o - 1];
  for (var s = new Array(4), o = 0; o < 4; o++)
    s[o] = $t();
  for (var c = new U("", e * 1.05 + 100 >> 0, e * 1.05 + 100 >> 0), o = e - 1; o >= 0; o--)
    s[o % 4] = Be(s[o % 4], c, i[n[o]], r[n[o]], 12);
  for (var o = 3; o >= 0; o--)
    zt(s[o], c);
  return _.Buffer.concat(
    [
      t.buf.slice(0, t.pos),
      c.buf.slice(c.pos, c.length)
    ],
    t.pos + c.length - c.pos
  );
}
function pn(n, e, t, r) {
  for (var a = 0; a < 256; a++) {
    e[a] = new Array(256), t[a] = new Array(256);
    for (var i = 0; i < 256; i++)
      e[a][i] = 0;
  }
  for (var o = Tt(n), a = 0; a < 256; a++)
    if (o[a]) {
      for (var s = 0, i = 0; i < 256; i++)
        o[i] && (s > 0 ? s-- : (e[a][i] = n.ReadUint7(), e[a][i] == 0 && (s = n.ReadByte())));
      it(e[a], r), t[a][0] = 0;
      for (var i = 0; i < 256; i++)
        t[a][i + 1] = t[a][i] + e[a][i];
    }
}
function _n(n, e) {
  var t = n.ReadByte(), r = t >> 4, a = n;
  if (t & 1)
    var i = n.ReadUint7(), o = n.ReadUint7(), t = new U(n.ReadData(o)), a = new U(at(t, i));
  var s = new Array(256), c = new Array(256);
  pn(a, s, c, r);
  for (var l = new Array(256), f = 0; f < 256; f++)
    l[f] = Dt(c[f], r);
  for (var d = new Array(4), h = new Array(4), u = 0; u < 4; u++)
    d[u] = n.ReadUint32(), h[u] = 0;
  for (var g = new _.Buffer.allocUnsafe(e), m = Math.floor(e / 4), f = 0; f < m; f++)
    for (var u = 0; u < 4; u++) {
      var v = Te(d[u], r), w = l[h[u]][v];
      g[f + u * m] = w, d[u] = He(d[u], c[h[u]][w], s[h[u]][w], r), d[u] = Oe(n, d[u]), h[u] = w;
    }
  for (f = 4 * f; f < e; ) {
    var v = Te(d[3], r), w = an(c[h[3]], v);
    g[f++] = w, d[3] = He(d[3], c[h[3]][w], s[h[3]][w], r), d[3] = Oe(n, d[3]), h[3] = w;
  }
  return g;
}
function Sn(n, e, t) {
  for (var r = 0; r < 256; r++) {
    t[r] = 0;
    for (var a = 0; a < 256; a++)
      e[r][a] = 0;
  }
  for (var i = 0, r = 0; r < n.length; r++)
    t[i]++, e[i][n[r]]++, i = n[r];
  t[i]++, e[0][n[1 * (n.length >> 2)]]++, e[0][n[2 * (n.length >> 2)]]++, e[0][n[3 * (n.length >> 2)]]++, t[0] += 3;
}
function An(n, e, t) {
  for (var r = 0; r < 256; r++)
    if (e[r]) {
      var a = Math.ceil(Math.log2(e[r]));
      a > t && (a = t), Qe(n[r], a);
    }
}
function bn(n, e, t) {
  for (var r = 0; r < 256; r++)
    e[r] && it(n[r], t);
}
function Cn(n, e, t) {
  Ht(n, t);
  for (var r = 0; r < 256; r++)
    if (t[r]) {
      for (var a = 0, i = 0; i < 256; i++)
        if (t[i]) {
          if (a)
            a--;
          else if (n.WriteUint7(e[r][i]), !e[r][i]) {
            for (var o = i + 1; o < 256; o++)
              if (t[o])
                if (e[r][o] == 0)
                  a++;
                else
                  break;
            n.WriteByte(a);
          }
        }
    }
}
function Bn(n) {
  const e = n.length;
  for (var t = new U("", 0, 257 * 257 * 3 + 9), r = new Array(256), a = new Array(256), i = new Array(256), o = 0; o < 256; o++)
    a[o] = new Array(256), i[o] = new Array(256);
  var s = 12;
  Sn(n, a, r), An(a, r, s);
  var c = new U("", 0, 257 * 257 * 3 + 9);
  Cn(c, a, r);
  var l = ot(c.buf.slice(0, c.pos));
  l.length < c.pos ? (t.WriteByte(1 | s << 4), t.WriteUint7(c.pos), t.WriteUint7(l.length), t.WriteData(l, l.length)) : (t.WriteByte(0 | s << 4), t.WriteData(c.buf, c.pos)), bn(a, r, s);
  for (var o = 0; o < 256; o++)
    if (r[o]) {
      i[o][0] = 0;
      for (var f = 1; f < 256; f++)
        i[o][f] = i[o][f - 1] + a[o][f - 1];
    }
  for (var d = new Array(4), h = new Array(4), f = 0; f < 4; f++)
    d[f] = $t(), h[f] = 0;
  for (var u = new U("", e * 1.05 + 100 >> 0, e * 1.05 + 100 >> 0), g = Math.floor(e / 4), m = new Array(4), v = new Array(4), f = 0; f < 4; f++)
    m[f] = (f + 1) * g - 2, v[f] = n[m[f] + 1];
  v[3] = n[e - 1];
  for (var o = e - 2; o > 4 * g - 2; o--)
    d[3] = Be(d[3], u, i[n[o]][v[3]], a[n[o]][v[3]], s), v[3] = n[o];
  for (; m[0] >= 0; )
    for (var f = 3; f >= 0; f--) {
      var w = n[m[f]];
      d[f] = Be(d[f], u, i[w][v[f]], a[w][v[f]], s), v[f] = w, m[f]--;
    }
  for (var f = 3; f >= 0; f--)
    d[f] = Be(d[f], u, i[0][v[f]], a[0][v[f]], s);
  for (var o = 3; o >= 0; o--)
    zt(d[o], u);
  return _.Buffer.concat(
    [
      t.buf.slice(0, t.pos),
      u.buf.slice(u.pos, u.length)
    ],
    t.pos + u.length - u.pos
  );
}
var Ot = { decode: vn, encode: Ge }, Gt = class {
  constructor(e) {
    this.low = 0, this.range = 4294967295, this.code = 0, this.FFnum = 0, this.carry = 0, this.cache = 0;
  }
  RangeStartDecode(e) {
    for (var t = 0; t < 5; t++)
      this.code = (this.code << 8) + e.ReadByte();
    this.code &= 4294967295, this.code >>>= 0;
  }
  RangeGetFrequency(e) {
    return this.range = Math.floor(this.range / e), Math.floor(this.code / this.range);
  }
  RangeDecode(e, t, r, a) {
    for (this.code -= t * this.range, this.range *= r; this.range < 1 << 24; )
      this.range *= 256, this.code = this.code * 256 + e.ReadByte();
  }
  RangeShiftLow(e) {
    if (this.low < 4278190080 | this.carry) {
      for (e.WriteByte(this.cache + this.carry); this.FFnum; )
        e.WriteByte(this.carry - 1), this.FFnum--;
      this.cache = this.low >>> 24, this.carry = 0;
    } else
      this.FFnum++;
    this.low <<= 8, this.low >>>= 0;
  }
  RangeEncode(e, t, r, a) {
    var i = this.low;
    for (this.range = Math.floor(this.range / a), this.low += t * this.range, this.low >>>= 0, this.range *= r, this.low < i && (this.carry != 0 && console.log("ERROR: Multiple carry"), this.carry = 1); this.range < 1 << 24; )
      this.range *= 256, this.RangeShiftLow(e);
  }
  RangeFinishEncode(e) {
    for (var t = 0; t < 5; t++)
      this.RangeShiftLow(e);
  }
};
const vt = 65536 - 17, ye = 16;
var Qt = class {
  constructor(e = 256) {
    this.total_freq = e, this.max_sym = e - 1, this.S = new Array(), this.F = new Array();
    for (var t = 0; t <= this.max_sym; t++)
      this.S[t] = t, this.F[t] = 1;
  }
  ModelDecode(e, t) {
    for (var r = t.RangeGetFrequency(this.total_freq), a = 0, i = 0; a + this.F[i] <= r; )
      a += this.F[i++];
    t.RangeDecode(e, a, this.F[i], this.total_freq), this.F[i] += ye, this.total_freq += ye, this.total_freq > vt && this.ModelRenormalise();
    var o = this.S[i];
    if (i > 0 && this.F[i] > this.F[i - 1]) {
      var s = this.F[i];
      this.F[i] = this.F[i - 1], this.F[i - 1] = s, s = this.S[i], this.S[i] = this.S[i - 1], this.S[i - 1] = s;
    }
    return o;
  }
  ModelRenormalise() {
    this.total_freq = 0;
    for (var e = 0; e <= this.max_sym; e++)
      this.F[e] -= Math.floor(this.F[e] / 2), this.total_freq += this.F[e];
  }
  ModelEncode(e, t, o) {
    for (var a = 0, i = 0; this.S[i] != o; i++)
      a += this.F[i];
    t.RangeEncode(e, a, this.F[i], this.total_freq), this.F[i] += ye, this.total_freq += ye, this.total_freq > vt && this.ModelRenormalise();
    var o = this.S[i];
    if (i > 0 && this.F[i] > this.F[i - 1]) {
      var s = this.F[i];
      this.F[i] = this.F[i - 1], this.F[i - 1] = s, s = this.S[i], this.S[i] = this.S[i - 1], this.S[i - 1] = s;
    }
  }
}, K = {};
K.array = function(n) {
  var e = 0, t = 0, r = [0, 1, 3, 7, 15, 31, 63, 127, 255];
  return function(a) {
    for (var i = 0; a > 0; ) {
      var o = 8 - e;
      a >= o ? (i <<= o, i |= r[o] & n[t++], e = 0, a -= o) : (i <<= a, i |= (n[t] & r[a] << 8 - a - e) >> 8 - a - e, e += a, a = 0);
    }
    return i;
  };
};
K.simple = function(n) {
  var e = K.header(n), t, r, a = [], i = 0;
  do
    r = K.decompress(n, e), r != -1 && (a.push(r), i += r.byteLength);
  while (r != -1);
  t = new Uint8Array(i), i = 0;
  for (var o = 0; o < a.length; ++o)
    r = a[o], t.set(r, i), i += r.byteLength;
  return t;
};
K.header = function(n) {
  if (n(8 * 3) != 4348520)
    throw "No magic number found";
  var e = n(8) - 48;
  if (e < 1 || e > 9)
    throw "Not a BZIP archive";
  return e;
};
K.decompress = function(n, e, t) {
  for (var r = 20, a = 258, i = 0, o = 1, s = 50, c = 1e5 * 9, l = "", f = 0; f < 6; f++)
    l += n(8).toString(16);
  if (l == "177245385090")
    return -1;
  if (l != "314159265359")
    throw "eek not valid bzip data";
  if (n(32), n(1))
    throw "unsupported obsolete version";
  var d = n(24);
  if (d > c)
    throw "Initial position larger than buffer size";
  var h = n(16), u = new Uint8Array(256), g = 0;
  for (f = 0; f < 16; f++)
    if (h & 1 << 15 - f) {
      var m = n(16);
      for (q = 0; q < 16; q++)
        m & 1 << 15 - q && (u[g++] = 16 * f + q);
    }
  var v = n(3);
  if (v < 2 || v > 6)
    throw "another error";
  var w = n(15);
  if (w == 0)
    throw "meh";
  for (var A = [], f = 0; f < v; f++)
    A[f] = f;
  for (var M = new Uint8Array(32768), f = 0; f < w; f++) {
    for (var q = 0; n(1); q++)
      if (q >= v)
        throw "whoops another error";
    var y = A[q];
    A.splice(q, 1), A.splice(0, 0, y), M[f] = y;
  }
  for (var T = g + 2, p = [], q = 0; q < v; q++) {
    var B = new Uint8Array(a), k = new Uint8Array(r + 1);
    h = n(5);
    for (var f = 0; f < T; f++) {
      for (; ; ) {
        if (h < 1 || h > r)
          throw "I gave up a while ago on writing error messages";
        if (!n(1))
          break;
        n(1) ? h-- : h++;
      }
      B[f] = h;
    }
    var I, L;
    I = L = B[0];
    for (var f = 1; f < T; f++)
      B[f] > L ? L = B[f] : B[f] < I && (I = B[f]);
    var R;
    R = p[q] = {}, R.permute = new Uint32Array(a), R.limit = new Uint32Array(r + 1), R.base = new Uint32Array(r + 1), R.minLen = I, R.maxLen = L;
    for (var W = R.base.subarray(1), F = R.limit.subarray(1), b = 0, f = I; f <= L; f++)
      for (var h = 0; h < T; h++)
        B[h] == f && (R.permute[b++] = h);
    for (f = I; f <= L; f++)
      k[f] = F[f] = 0;
    for (f = 0; f < T; f++)
      k[B[f]]++;
    for (b = h = 0, f = I; f < L; f++)
      b += k[f], F[f] = b - 1, b <<= 1, W[f + 1] = b - (h += k[f]);
    F[L] = b + k[L] - 1, W[I] = 0;
  }
  for (var D = new Uint32Array(256), f = 0; f < 256; f++)
    A[f] = f;
  var $, N, T, Me;
  $ = N = T = Me = 0;
  for (var J = new Uint32Array(c); ; ) {
    if (!T--) {
      if (T = s - 1, Me >= w)
        throw "meow i'm a kitty, that's an error";
      R = p[M[Me++]], W = R.base.subarray(1), F = R.limit.subarray(1);
    }
    for (f = R.minLen, q = n(f); ; ) {
      if (f > R.maxLen)
        throw "rawr i'm a dinosaur";
      if (q <= F[f])
        break;
      f++, q = q << 1 | n(1);
    }
    if (q -= W[f], q < 0 || q >= a)
      throw "moo i'm a cow";
    var ne = R.permute[q];
    if (ne == i || ne == o) {
      $ || ($ = 1, h = 0), ne == i ? h += $ : h += 2 * $, $ <<= 1;
      continue;
    }
    if ($) {
      if ($ = 0, N + h >= c)
        throw "Boom.";
      for (y = u[A[0]], D[y] += h; h--; )
        J[N++] = y;
    }
    if (ne > g)
      break;
    if (N >= c)
      throw "I can't think of anything. Error";
    f = ne - 1, y = A[f], A.splice(f, 1), A.splice(0, 0, y), y = u[y], D[y]++, J[N++] = y;
  }
  if (d < 0 || d >= N)
    throw "I'm a monkey and I'm throwing something at someone, namely you";
  for (var q = 0, f = 0; f < 256; f++)
    m = q + D[f], D[f] = q, q = m;
  for (var f = 0; f < N; f++)
    y = J[f] & 255, J[D[y]] |= f << 8, D[y]++;
  var Q = 0, V = 0, Fe = 0;
  N && (Q = J[d], V = Q & 255, Q >>= 8, Fe = -1), N = N;
  var qe = new Uint8Array(c), Pe, Le, ke, dt = 0;
  for (t || (t = 1 / 0); N; ) {
    for (N--, Le = V, Q = J[Q], V = Q & 255, Q >>= 8, Fe++ == 3 ? (Pe = V, ke = Le, V = -1) : (Pe = 1, ke = V); Pe--; )
      if (qe[dt++] = ke, !--t)
        return qe;
    V != Le && (Fe = 0);
  }
  return qe.subarray(0, dt);
};
var Rn = K;
const O = Gt, pe = ve, z = Qt, Ne = Rn, mt = 1, En = 4, wt = 8, yt = 16, xn = 32, pt = 64, _e = 128;
var Vt = class {
  decode(e) {
    return this.stream = new pe(e), this.decodeStream(this.stream);
  }
  decodeStream(e, t = 0) {
    var r = this.stream.ReadByte();
    r & yt || (t = this.stream.ReadUint7());
    var a = t, i = r & mt;
    if (r & wt)
      return this.decodeStripe(this.stream, t);
    if (r & _e) {
      var o;
      [o, a] = this.decodePackMeta(this.stream);
    }
    if (r & xn)
      var s = this.decodeCat(this.stream, a);
    else if (r & En)
      var s = this.decodeExt(this.stream, a);
    else if (r & pt)
      var s = i ? this.decodeRLE1(this.stream, a) : this.decodeRLE0(this.stream, a);
    else
      var s = i ? this.decode1(this.stream, a) : this.decode0(this.stream, a);
    return r & _e && (s = this.decodePack(s, o, t)), s;
  }
  encode(e, t) {
    if (this.stream = new pe("", 0, e.length * 1.1 + 100), this.stream.WriteByte(t), t & yt || this.stream.WriteUint7(e.length), t & wt)
      return _.Buffer.concat([
        this.stream.buf.slice(0, this.stream.pos),
        this.encodeStripe(this.stream, e, t >> 8)
      ]);
    var r = t & mt, a = e.length, i;
    return t & _e && ([i, e, a] = this.encodePack(e)), t & _e && this.stream.WriteStream(i), t & pt ? r ? this.encodeRLE1(e, a, this.stream) : this.encodeRLE0(e, a, this.stream) : r ? this.encode1(e, a, this.stream) : this.encode0(e, a, this.stream);
  }
  //----------------------------------------------------------------------
  // Order-0 codec
  decode0(e, t) {
    var r = new _.Buffer.allocUnsafe(t), a = e.ReadByte();
    a == 0 && (a = 256);
    var i = new z(a), o = new O(e);
    o.RangeStartDecode(e);
    for (var s = 0; s < t; s++)
      r[s] = i.ModelDecode(e, o);
    return r;
  }
  encode0(e, t, r) {
    for (var a = 0, i = 0; i < t; i++)
      a < e[i] && (a = e[i]);
    a++;
    var o = new z(a);
    r.WriteByte(a);
    for (var s = new O(r), i = 0; i < t; i++)
      o.ModelEncode(r, s, e[i]);
    return s.RangeFinishEncode(r), r.buf.slice(0, r.pos);
  }
  //----------------------------------------------------------------------
  // Order-1 codec
  decode1(e, t) {
    var r = new _.Buffer.allocUnsafe(t), a = e.ReadByte();
    a == 0 && (a = 256);
    for (var i = new Array(a), o = 0; o < a; o++)
      i[o] = new z(a);
    var s = new O(e);
    s.RangeStartDecode(e);
    for (var c = 0, o = 0; o < t; o++)
      r[o] = i[c].ModelDecode(e, s), c = r[o];
    return r;
  }
  encode1(e, t, r) {
    for (var a = 0, i = 0; i < t; i++)
      a < e[i] && (a = e[i]);
    a++;
    for (var o = new Array(a), i = 0; i < a; i++)
      o[i] = new z(a);
    r.WriteByte(a);
    for (var s = new O(r), c = 0, i = 0; i < t; i++)
      o[c].ModelEncode(r, s, e[i]), c = e[i];
    return s.RangeFinishEncode(r), r.buf.slice(0, r.pos);
  }
  //----------------------------------------------------------------------
  // External codec
  decodeExt(e, t) {
    var r = new _.Buffer.allocUnsafe(t), a = Ne.array(e.buf.slice(e.pos)), i = Ne.header(a), o = 0;
    do {
      var s = Ne.decompress(a, i);
      s != -1 && (_.Buffer.from(s).copy(r, o), o += s.length, i -= s.length);
    } while (s != -1);
    return r;
  }
  encodeExt(e, t) {
  }
  //----------------------------------------------------------------------
  // Order-0 RLE codec
  decodeRLE0(e, t) {
    var r = new _.Buffer.allocUnsafe(t), a = e.ReadByte();
    a == 0 && (a = 256);
    for (var i = new z(a), o = new Array(258), c = 0; c <= 257; c++)
      o[c] = new z(4);
    var s = new O(e);
    s.RangeStartDecode(e);
    for (var c = 0; c < t; ) {
      r[c] = i.ModelDecode(e, s);
      for (var l = o[r[c]].ModelDecode(e, s), f = l, d = 256; l == 3; )
        l = o[d].ModelDecode(e, s), d = 257, f += l;
      for (var h = 1; h <= f; h++)
        r[c + h] = r[c];
      c += f + 1;
    }
    return r;
  }
  encodeRLE0(e, t, r) {
    for (var a = 0, c = 0; c < t; c++)
      a < e[c] && (a = e[c]);
    a++;
    for (var i = new z(a), o = new Array(258), c = 0; c <= 257; c++)
      o[c] = new z(4);
    r.WriteByte(a);
    for (var s = new O(r), c = 0; c < t; ) {
      i.ModelEncode(r, s, e[c]);
      for (var l = 1; c + l < t && e[c + l] == e[c]; )
        l++;
      l--;
      var f = e[c];
      e[c], c += l + 1;
      var d = l >= 3 ? 3 : l;
      for (o[f].ModelEncode(r, s, d), l -= d, f = 256; d == 3; )
        d = l >= 3 ? 3 : l, o[f].ModelEncode(r, s, d), f = 257, l -= d;
    }
    return s.RangeFinishEncode(r), r.buf.slice(0, r.pos);
  }
  //----------------------------------------------------------------------
  // Order-1 RLE codec
  decodeRLE1(e, t) {
    var r = new _.Buffer.allocUnsafe(t), a = e.ReadByte();
    a == 0 && (a = 256);
    for (var i = new Array(a), l = 0; l < a; l++)
      i[l] = new z(a);
    for (var o = new Array(258), l = 0; l <= 257; l++)
      o[l] = new z(4);
    var s = new O(e);
    s.RangeStartDecode(e);
    for (var c = 0, l = 0; l < t; ) {
      r[l] = i[c].ModelDecode(e, s), c = r[l];
      for (var f = o[r[l]].ModelDecode(e, s), d = f, h = 256; f == 3; )
        f = o[h].ModelDecode(e, s), h = 257, d += f;
      for (var u = 1; u <= d; u++)
        r[l + u] = r[l];
      l += d + 1;
    }
    return r;
  }
  encodeRLE1(e, t, r) {
    for (var a = 0, c = 0; c < t; c++)
      a < e[c] && (a = e[c]);
    a++;
    for (var i = new Array(a), c = 0; c < a; c++)
      i[c] = new z(a);
    for (var o = new Array(258), c = 0; c <= 257; c++)
      o[c] = new z(4);
    r.WriteByte(a);
    for (var s = new O(r), c = 0, l = 0; c < t; ) {
      i[l].ModelEncode(r, s, e[c]);
      for (var f = 1; c + f < t && e[c + f] == e[c]; )
        f++;
      f--;
      var d = e[c];
      l = e[c], c += f + 1;
      var h = f >= 3 ? 3 : f;
      for (o[d].ModelEncode(r, s, h), f -= h, d = 256; h == 3; )
        h = f >= 3 ? 3 : f, o[d].ModelEncode(r, s, h), d = 257, f -= h;
    }
    return s.RangeFinishEncode(r), r.buf.slice(0, r.pos);
  }
  //----------------------------------------------------------------------
  // Pack method
  decodePackMeta(e) {
    this.nsym = e.ReadByte();
    for (var t = new Array(this.nsym), r = 0; r < this.nsym; r++)
      t[r] = e.ReadByte();
    var a = e.ReadUint7();
    return [t, a];
  }
  decodePack(e, t, r) {
    var a = new _.Buffer.allocUnsafe(r);
    if (this.nsym <= 1)
      for (var i = 0; i < r; i++)
        a[i] = t[0];
    else if (this.nsym <= 2)
      for (var i = 0, o = 0; i < r; i++) {
        if (i % 8 == 0)
          var s = e[o++];
        a[i] = t[s & 1], s >>= 1;
      }
    else if (this.nsym <= 4)
      for (var i = 0, o = 0; i < r; i++) {
        if (i % 4 == 0)
          var s = e[o++];
        a[i] = t[s & 3], s >>= 2;
      }
    else if (this.nsym <= 16)
      for (var i = 0, o = 0; i < r; i++) {
        if (i % 2 == 0)
          var s = e[o++];
        a[i] = t[s & 15], s >>= 4;
      }
    else
      return e;
    return a;
  }
  // Compute M array and return meta-data stream
  packMeta(e) {
    for (var t = new pe("", 0, 1024), r = new Array(256), a = 0; a < e.length; a++)
      r[e[a]] = 1;
    for (var i = 0, a = 0; a < 256; a++)
      r[a] && (r[a] = ++i);
    t.WriteByte(i);
    for (var a = 0; a < 256; a++)
      r[a] && (t.WriteByte(a), r[a]--);
    return [t, r, i];
  }
  encodePack(e) {
    var t, r, a;
    [t, r, a] = this.packMeta(e);
    var i = e.length, o = 0;
    if (a <= 1)
      return t.WriteUint7(0), [t, new _.Buffer.allocUnsafe(0), 0];
    if (a <= 2) {
      for (var s = new _.Buffer.allocUnsafe(Math.floor((i + 7) / 8)), o = 0, c = 0; o < (i & -8); o += 8, c++)
        s[c] = (r[e[o + 0]] << 0) + (r[e[o + 1]] << 1) + (r[e[o + 2]] << 2) + (r[e[o + 3]] << 3) + (r[e[o + 4]] << 4) + (r[e[o + 5]] << 5) + (r[e[o + 6]] << 6) + (r[e[o + 7]] << 7);
      if (o < i) {
        s[c] = 0;
        for (var l = 0; o < i; )
          s[c] |= r[e[o++]] << l, l++;
        c++;
      }
      return t.WriteUint7(c), [t, s, s.length];
    }
    if (a <= 4) {
      for (var s = new _.Buffer.allocUnsafe(Math.floor((i + 3) / 4)), o = 0, c = 0; o < (i & -4); o += 4, c++)
        s[c] = (r[e[o + 0]] << 0) + (r[e[o + 1]] << 2) + (r[e[o + 2]] << 4) + (r[e[o + 3]] << 6);
      if (o < i) {
        s[c] = 0;
        for (var l = 0; o < i; )
          s[c] |= r[e[o++]] << l, l += 2;
        c++;
      }
      return t.WriteUint7(c), [t, s, s.length];
    }
    if (a <= 16) {
      for (var s = new _.Buffer.allocUnsafe(Math.floor((i + 1) / 2)), o = 0, c = 0; o < (i & -2); o += 2, c++)
        s[c] = (r[e[o + 0]] << 0) + (r[e[o + 1]] << 4);
      return o < i && (s[c++] = r[e[o++]]), t.WriteUint7(c), [t, s, s.length];
    }
    return t.WriteUint7(e.length), [t, e, e.length];
  }
  //----------------------------------------------------------------------
  // STRIPE method
  encodeStripe(e, t, r) {
    r == 0 && (r = 4);
    for (var a = new Array(r), i = new Array(r), o = 0; o < r; o++)
      i[o] = Math.floor(t.length / r) + (t.length % r > o), a[o] = new Array(i[o]);
    for (var s = 0, c = 0; c < t.length; c += r, s++)
      for (var l = 0; l < r; l++)
        s < a[l].length && (a[l][s] = t[c + l]);
    for (var f = new Array(r), d = 0, o = 0; o < r; o++) {
      var h = this.encode(a[o], 0), u = this.encode(a[o], 1);
      f[o] = u.length < h.length ? u : h, d += f[o].length;
    }
    var g = new pe("", 0, d + 5 * r + 1);
    g.WriteByte(r);
    for (var o = 0; o < r; o++)
      g.WriteUint7(f[o].length);
    for (var o = 0; o < r; o++)
      g.WriteData(f[o], f[o].length);
    return g.buf.slice(0, g.buf.pos);
  }
  decodeStripe(e, t) {
    for (var r = e.ReadByte(), a = new Array(r), i = new Array(r), o = 0; o < r; o++)
      a[o] = e.ReadUint7();
    for (var s = new Array(r), o = 0; o < r; o++)
      i[o] = Math.floor(t / r) + (t % r > o), s[o] = this.decodeStream(e, i[o]);
    for (var c = new _.Buffer.allocUnsafe(t), o = 0; o < r; o++)
      for (var l = 0; l < i[o]; l++)
        c[l * r + o] = s[o][l];
    return c;
  }
  //----------------------------------------------------------------------
  // Cat method
  decodeCat(e, t) {
    for (var r = new _.Buffer.allocUnsafe(t), a = 0; a < t; a++)
      r[a] = e.ReadByte();
    return r;
  }
};
const Kt = ve, H = Qt, Xt = Gt;
function Re(n, e, t) {
  for (var r = 0, a = 0, i = -1, o = new Array(1024); a < t; ) {
    var s = n.ReadByte();
    if (o[r++] = s, a += s, s == i) {
      var c = n.ReadByte();
      for (a += s * c; c--; )
        o[r++] = s;
    }
    i = s;
  }
  var l = 0;
  for (r = 0, a = 0; a < t; ) {
    var f = 0;
    do {
      var d = o[r++];
      f += d;
    } while (d == 255);
    for (; f--; )
      e[a++] = l;
    l++;
  }
}
const st = 2, Yt = 4, Zt = 8, Ve = 16, Ke = 32, Xe = 64, Ye = 128, Ze = 1, Je = 2, In = 4;
function Mn(n, e, t) {
  var r = n.context;
  return e.qctx = (e.qctx << n.qshift) + n.qtab[t], r += (e.qctx & (1 << n.qbits) - 1) << n.qloc, n.do_pos && (r += n.ptab[Math.min(e.p, 1023)] << n.ploc), n.do_delta && (r += n.dtab[Math.min(e.delta, 255)] << n.dloc, e.delta += e.prevq != t ? 1 : 0, e.prevq = t), n.do_sel && (r += e.s << n.sloc), e.p--, r & 65535;
}
function Fn(n) {
  var e = {};
  e.context = n.ReadUint16(), e.pflags = n.ReadByte(), e.do_dedup = e.pflags & st, e.fixed_len = e.pflags & Yt, e.do_sel = e.pflags & Zt, e.do_qmap = e.pflags & Ve, e.do_pos = e.pflags & Ke, e.do_delta = e.pflags & Xe, e.do_qtab = e.pflags & Ye, e.max_sym = n.ReadByte();
  var t = n.ReadByte();
  if (e.qbits = t >> 4, e.qshift = t & 15, t = n.ReadByte(), e.qloc = t >> 4, e.sloc = t & 15, t = n.ReadByte(), e.ploc = t >> 4, e.dloc = t & 15, e.qmap = new Array(256), e.pflags & Ve)
    for (var r = 0; r < e.max_sym; r++)
      e.qmap[r] = n.ReadByte();
  else
    for (var r = 0; r < 256; r++)
      e.qmap[r] = r;
  if (e.qtab = new Array(1024), e.qbits > 0 && e.pflags & Ye)
    Re(n, e.qtab, 256);
  else
    for (var r = 0; r < 256; r++)
      e.qtab[r] = r;
  return e.ptab = new Array(1024), e.pflags & Ke && Re(n, e.ptab, 1024), e.dtab = new Array(256), e.pflags & Xe && Re(n, e.dtab, 256), e;
}
function qn(n) {
  var e = {
    max_sym: 0
  }, t = n.ReadByte();
  if (t != 5) {
    console.error("Invalid FQZComp version number");
    return;
  }
  var r = n.ReadByte(), a = r & Ze ? n.ReadByte() : 1, i = r.nparam > 1 ? r.nparam - 1 : 0, o = new Array(256);
  if (r & Je)
    i = n.ReadByte(), Re(n, o, 256);
  else {
    for (var s = 0; s < a; s++)
      o[s] = s;
    for (; s < 256; s++)
      o[s] = a - 1;
  }
  e.do_rev = r & In, e.stab = o, e.max_sel = i, e.params = new Array(e.nparam);
  for (var c = 0; c < a; c++)
    e.params[c] = Fn(n), e.max_sym < e.params[c].max_sym && (e.max_sym = e.params[c].max_sym);
  return e;
}
function Pn(n) {
  var e = {};
  e.qual = new Array(65536);
  for (var t = 0; t < 65536; t++)
    e.qual[t] = new H(n.max_sym + 1);
  e.len = new Array(4);
  for (var t = 0; t < 4; t++)
    e.len[t] = new H(256);
  return e.rev = new H(2), e.dup = new H(2), n.max_sel > 0 && (e.sel = new H(n.max_sel + 1)), e;
}
function Ln(n, e, t, r, a, i) {
  t.max_sel > 0 ? a.s = r.sel.ModelDecode(n, e) : a.s = 0, a.x = t.stab[a.s];
  var o = t.params[a.x];
  if (o.fixed_len >= 0) {
    var s = r.len[0].ModelDecode(n, e);
    s |= r.len[1].ModelDecode(n, e) << 8, s |= r.len[2].ModelDecode(n, e) << 16, s |= r.len[3].ModelDecode(n, e) << 24, o.fixed_len > 0 && (o.fixed_len = -s);
  } else
    s = -o.fixed_len;
  a.len = s, t.do_rev && (i[a.rec] = r.rev.ModelDecode(n, e)), a.is_dup = 0, o.pflags & st && r.dup.ModelDecode(n, e) && (a.is_dup = 1), a.p = s, a.delta = 0, a.qctx = 0, a.prevq = 0, a.rec++;
}
function kn(n, e) {
  var t = n.ReadUint7(), r = qn(n);
  if (r) {
    var a = r.params, i = new Array(e.length), o = Pn(r), s = new Xt(n);
    s.RangeStartDecode(n);
    for (var c = new _.Buffer.allocUnsafe(t), l = {
      qctx: 0,
      // Qual-only sub-context
      prevq: 0,
      // Previous quality value
      delta: 0,
      // Running delta (q vs prevq)
      p: 0,
      // Number of bases left in current record
      s: 0,
      // Current parameter selector value (0 if unused)
      x: 0,
      // "stab" tabulated copy of s
      len: 0,
      // Length of current string
      is_dup: 0,
      // This string is a duplicate of last
      rec: 0
      // Record number
    }, f = 0; f < t; ) {
      if (l.p == 0) {
        if (Ln(n, s, r, o, l, i), l.is_dup > 0 && o.dup.ModelDecode(n, s)) {
          for (var d = 0; d < len; d++)
            c[f + d] = c[f + d - l.len];
          f += l.len, l.p = 0;
          continue;
        }
        e.push(l.len);
        var a = r.params[l.x], h = a.context;
      }
      var u = o.qual[h].ModelDecode(n, s);
      c[f++] = a.qmap[u], h = Mn(a, l, u);
    }
    return r.do_rev && Un(c, t, i, e), c;
  }
}
function Un(n, e, t, r) {
  for (var a = 0, i = 0; i < e; ) {
    if (t[a])
      for (var o = 0, s = r[a] - 1; o < s; ) {
        var c = n[i + o];
        n[i + o] = n[i + s], n[i + s] = c, o++, s--;
      }
    i += r[a++];
  }
}
function Nn(n, e) {
  var t = new Kt(n);
  return kn(t, e);
}
function Dn(n, e, t, r) {
  for (var a = t[0], i = 0; i < t.length && t[i] == a; i++)
    ;
  for (var o = i == t.length ? 1 : 0, s = 0, c = 0, i = 0; i < 256; i++)
    r[0][i] = 0;
  for (var l = 0, f = 0, i = 0; i < n.length; i++)
    f == 0 && (f = e[l < e.length - 1 ? l++ : l]), r[0][n[i]]++, f--;
  for (var i = 0; i < 256; i++)
    r[0][i] && (c < i && (c = i), s++);
  var d = 5, h = 0;
  return s <= 16 && (h = 1, s <= 2 ? d = 1 : s <= 4 ? d = 2 : s <= 8 ? d = 3 : d = 4), [{
    qbits: 8 + (d > 4),
    qshift: d,
    qloc: 7,
    pbits: 7,
    pshift: e[0] > 128 ? 1 : 0,
    ploc: 0,
    dbits: d > 4 ? 0 : 1,
    dshift: 3,
    dloc: 15,
    // NB: Also useful as a way of embedding sel and doing sel
    // specific contexts. Identical bar context. Eg 0<<15 or 1<<15.
    sbits: 0,
    sloc: 15,
    do_stab: 0,
    context: 0,
    max_sym: c,
    nsym: s,
    do_qmap: h,
    do_dedup: 0,
    fixed_len: e.length == 1 ? 1 : 0,
    do_sel: 0,
    do_rev: 0,
    do_pos: 1,
    do_delta: d <= 4 ? 1 : 0,
    do_qtab: 0,
    // Override above with some attempt at using selectors
    // when the q_dirs are specific and non-fixed.
    qbits: 8 + (d > 4) - (o == 0),
    sbits: 1,
    sloc: 15 - (d <= 4),
    // read1 vs read2
    do_stab: 1,
    do_sel: 1
    //	     // q4+dir: 7245769 with, 7353962 without. 1.5% saving
    //	     qbits:     6,
    //	     dbits:     2,
    //	     dshift:    2,
    //	     dloc:      13,
    //	     sbits:     1,
    //	     sloc:      15,
    //	     do_stab:   1,
    //	     do_sel:    1,
    // with 20 bits of context, q40 = 31741545
    // qbits 10, dbits 2, pbits 7, sbits 1
  }];
}
function Se(n, e, t) {
  for (var r = 0, a = 0, i = new Array(t * 2), o = 0; r < t; ) {
    for (var s = r; r < t && e[r] == a; )
      r++;
    var c = r - s;
    do {
      var l = Math.min(255, c);
      i[o++] = l, c -= l;
    } while (l == 255);
    a++;
  }
  var f = -1, d = new Array(t * 2), h = 0;
  for (r = 0; r < o; ) {
    var u = i[r++];
    if (d[h++] = u, u == f) {
      for (var s = r; r < o && i[r] == f && r - s < 255; )
        r++;
      d[h++] = r - s;
    } else
      f = u;
  }
  n.WriteData(d, h);
}
function $n(n, e, t, r, a, i, o) {
  for (var s = [
    0,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7,
    7
  ], c = 0; c < e.length; c++)
    o[c] = c;
  for (; c < 256; c++)
    o[c] = e.length - 1;
  n.WriteByte(5);
  var l = (e.length > 1 ? Ze : 0) | (e[0].do_stab ? Je : 0);
  if (n.WriteByte(l), l & Ze && n.WriteByte(e.length), l & Je) {
    var f = 1 << e[0].sbits;
    f > 0 && f--, n.WriteByte(f), Se(n, o, 256);
  }
  for (var d = 0; d < e.length; d++) {
    if (n.WriteUint16(e[d].context), n.WriteByte((e[d].do_qtab ? Ye : 0) | // FLAG
    (e[d].do_delta ? Xe : 0) | (e[d].do_pos ? Ke : 0) | (e[d].do_qmap ? Ve : 0) | (e[d].do_sel ? Zt : 0) | (e[d].fixed_len ? Yt : 0) | (e[d].do_dedup ? st : 0)), e[d].do_qmap ? n.WriteByte(e[d].nsym) : n.WriteByte(e[d].max_sym), n.WriteByte(e[d].qbits << 4 | e[d].qshift), n.WriteByte(e[d].qloc << 4 | e[d].sloc), n.WriteByte(e[d].ploc << 4 | e[d].dloc), e[d].do_qmap) {
      e[d].max_sym = e[d].nsym;
      for (var h = 0, c = 0; c < 256; c++)
        t[d][c] && (n.WriteByte(c), t[d][c] = h++);
      for (; h < e[d].nsym; h++)
        n.WriteByte(0);
    } else
      for (var c = 0; c < 256; c++)
        t[d][c] = c;
    if (e[d].qbits > 0) {
      for (var c = 0; c < 256; c++)
        r[d][c] = c;
      e[d].do_qtab && Se(n, r[d], 256);
    }
    if (e[d].pbits > 0) {
      for (var c = 0; c < 1024; c++)
        a[d][c] = Math.min((1 << e[d].pbits) - 1, c >> e[d].pshift);
      Se(n, a[d], 1024);
    }
    if (e[d].dbits > 0) {
      for (var c = 0; c < 256; c++)
        s[c] > (1 << e[d].dbits) - 1 && (s[c] = (1 << e[d].dbits) - 1);
      for (var c = 0; c < 256; c++)
        i[d][c] = s[Math.min(s.length - 1, c >> e[d].dshift)];
      Se(n, i[d], 256);
    }
  }
  return n;
}
function zn(n, e, t, r, a, i, o, s, c, l) {
  var f = 1 << a[0].sbits;
  f > 0 && f--;
  for (var d = e.length, h = 0, w = 0; w < a.length; w++)
    h < a[w].max_sym && (h = a[w].max_sym);
  for (var u = new Array(65536), A = 0; A < 65536; A++)
    u[A] = new H(h + 1);
  for (var g = new Array(4), A = 0; A < 4; A++)
    g[A] = new H(256);
  new H(2), new H(2);
  for (var m = new H(f + 1), v = new Xt(e), w = 0, A = 0, M = 0; A < d; ) {
    if (w == 0) {
      var y = r[M];
      a[0].sbits > 0 && m.ModelEncode(n, v, y);
      var p = l[y], B = t[Math.min(t.length - 1, M++)];
      a[p].fixed_len ? a[p].fixed_len > 0 && (g[0].ModelEncode(n, v, B & 255), g[1].ModelEncode(n, v, B >> 8 & 255), g[2].ModelEncode(n, v, B >> 16 & 255), g[3].ModelEncode(n, v, B >> 24 & 255), a[p].fixed_len = -1) : (g[0].ModelEncode(n, v, B & 255), g[1].ModelEncode(n, v, B >> 8 & 255), g[2].ModelEncode(n, v, B >> 16 & 255), g[3].ModelEncode(n, v, B >> 24 & 255)), a[p].do_dedup && cr.exit(1), w = B;
      var k = 0, I = a[p].context, L = 0, R = 0;
    }
    var W = e[A++], F = i[p][W];
    u[I].ModelEncode(n, v, F), L = (L << a[p].qshift) + o[p][F], I = a[p].context, I += (L & (1 << a[p].qbits) - 1) << a[p].qloc, a[p].pbits > 0 && (I += s[p][Math.min(w, 1023)] << a[p].ploc), a[p].dbits > 0 && (I += c[p][Math.min(k, 255)] << a[p].dloc, k += R != F ? 1 : 0, R = F), a[p].do_sel && (I += y << a[p].sloc), I = I & 65535, w--;
  }
  return v.RangeFinishEncode(n), n.buf.slice(0, n.pos);
}
function Wn(n, e, t) {
  for (var r = new Array(2), a = new Array(2), i = new Array(2), o = new Array(2), s = new Array(256), c = 0; c < 2; c++)
    r[c] = new Array(256), a[c] = new Array(256), i[c] = new Array(1024), o[c] = new Array(256);
  var f = new Kt("", 0, n.length * 1.1 + 100);
  f.WriteUint7(n.length);
  var l = Dn(n, e, t, r), f = $n(f, l, r, a, i, o, s);
  return zn(f, n, e, t, l, r, a, i, o, s);
}
var Tn = { decode: Nn, encode: Wn };
const ee = ve, Jt = Ot, Hn = Vt;
var jt = new Hn();
const ce = 0, le = 1, de = 2, te = 3, er = 4, se = 5, je = 6, re = 7, he = 8, ue = 9, ft = 10, ge = 12;
function On(n, e, t, r) {
  for (var a = -1, i = new Array(256); !n.EOF(); ) {
    var o = n.ReadByte(), s = o & 128, c = o & 64, l = o & 63;
    if (s && (a++, i[a] = new Array(13)), l != ce && s) {
      var f = new Array(r - 1).fill(ft);
      i[a][ce] = new ee(_.Buffer.from([l].concat(f)));
    }
    if (c) {
      var d = n.ReadByte(), h = n.ReadByte();
      i[a][l] = new ee(i[d][h].buf);
    } else {
      var u = n.ReadUint7(), g = n.ReadData(u);
      t ? i[a][l] = jt.decode(g) : i[a][l] = Jt.decode(g), i[a][l] = new ee(i[a][l]);
    }
  }
  return i;
}
function _t(n, e) {
  for (var t = n + ""; t.length < e; )
    t = "0" + t;
  return t;
}
function Gn(n, e, t, r) {
  var a = n[0][ce].ReadByte(), i = n[0][a].ReadUint32(), o = r - i;
  if (a == se)
    return e[r] = e[o], t[r] = t[o], e[r];
  var s = 1;
  e[r] = "", t[r] = new Array(256);
  do {
    switch (a = n[s][ce].ReadByte(), a) {
      case de:
        t[r][s] = n[s][de].ReadChar();
        break;
      case le:
        t[r][s] = n[s][le].ReadString();
        break;
      case re:
        t[r][s] = n[s][re].ReadUint32();
        break;
      case te:
        var c = n[s][te].ReadUint32(), l = n[s][er].ReadByte();
        t[r][s] = _t(c, l);
        break;
      case he:
        t[r][s] = (t[o][s] >> 0) + n[s][he].ReadByte();
        break;
      case ue:
        var c = (t[o][s] >> 0) + n[s][ue].ReadByte(), l = t[o][s].length;
        t[r][s] = _t(c, l);
        break;
      case ft:
        t[r][s] = t[o][s];
        break;
      default:
        t[r][s] = "";
        break;
    }
    e[r] += t[r][s++];
  } while (a != ge);
  return e[r];
}
function Qn(r, e, t) {
  var r = new ee(r);
  r.ReadUint32();
  var a = r.ReadUint32(), i = r.ReadByte(), o = On(r, e, i, a), s = new Array(a), c = new Array(a), l = "";
  typeof t > "u" && (t = `
`);
  for (var f = 0; f < a; f++)
    l += Gn(o, s, c, f) + t;
  return l;
}
function Vn(n, e) {
  var t = n.toString();
  t[t.length - 1] == `
` && (t = t.substring(0, t.length - 1));
  var r = t.split(`
`), a = new ee("", 0, t.length * 2 + 1e4);
  a.WriteUint32(t.length), a.WriteUint32(r.length), a.WriteByte(e);
  for (var i = new Array(r.length), o = {}, s = new Array(256).fill(0), c = 0, l = 0, f = 0; f < r.length; f++) {
    var [d, h] = Zn(i, o, s, r[f], f);
    c < d && (c = d), l < h && (l = h);
  }
  for (var u = 0; u < c; u++) {
    for (var g = new Array(ge + 1), m = 0; m <= ge; m++)
      g[m] = new ee("", 0, r.length * l);
    Kn(g, i, u, r), Xn(g, u, e, a);
  }
  return a.buf.slice(0, a.pos);
}
function Kn(n, e, t, r, a, i) {
  for (var o = 0; o < r.length; o++)
    if (!(t > 0 && e[o][0].type == se) && e[o][t])
      switch (n[ce].WriteByte(e[o][t].type), e[o][t].type) {
        case je:
          n[je].WriteUint32(e[o][t].val);
          break;
        case se:
          n[se].WriteUint32(e[o][t].val);
          break;
        case le:
          n[le].WriteString(e[o][t].val);
          break;
        case de:
          n[de].WriteChar(e[o][t].val);
          break;
        case re:
          n[re].WriteUint32(e[o][t].val);
          break;
        case te:
          n[te].WriteUint32(e[o][t].val), n[er].WriteByte(e[o][t].val.length);
          break;
        case he:
          n[e[o][t].type].WriteByte(e[o][t].val);
          break;
        case ue:
          n[e[o][t].type].WriteByte(e[o][t].val);
          break;
      }
}
function Xn(n, e, t, r) {
  for (var a = 0; a <= ge; a++)
    if (!(n[a].pos <= 0)) {
      r.WriteByte(a + (a == 0 ? 128 : 0)), n[a] = n[a].buf.slice(0, n[a].pos);
      var i = Yn(n[a], t);
      r.WriteUint7(i.length), r.WriteData(i, i.length);
    }
}
function Yn(n, e) {
  var t = 1073741824, r, a = [0, 1, 64, 65, 128, 129, 193 + 8];
  for (var i in a) {
    var o = a[i];
    if (!(o & 1 && n.length < 100) && !(o & 8 && n.length % 4 != 0)) {
      try {
        var s = e ? jt.encode(n, o) : Jt.encode(n, o);
      } catch {
        var s = 0;
      }
      s && t > s.length && (t = s.length, r = s);
    }
  }
  return r;
}
function Zn(n, e, t, r, a) {
  var i = 0, o = a - 1;
  n[a] = new Array(256), e[r] ? n[a][0] = {
    type: se,
    val: a - e[r]
  } : n[a][0] = {
    type: je,
    val: a == 0 ? 0 : 1
  }, e[r] = a;
  for (var s = r.match(/([a-zA-Z0-9]{1,9})|([^a-zA-Z0-9]+)/g), c = 0; c < s.length; c++) {
    var l = c + 1, f = le, d = s[c];
    if (s[c].match(/^0+[0-9]*$/g) ? f = te : s[c].match(/^[0-9]+$/g) ? f = re : s[c].length == 1 && (f = de), o >= 0 && n[o][l]) {
      if (n[o][l].str == s[c])
        f = ft, d = "";
      else if (n[o][l].type == re || n[o][l].type == he) {
        var h = d - n[o][l].str;
        t[l]++, h >= 0 && h < 256 && t[l] > a / 2 && (f = he, d = h);
      } else if ((n[o][l].type == te || n[o][l].type == ue) && n[o][l].str.length == d.length) {
        var h = d - n[o][l].str;
        t[l]++, h >= 0 && h < 256 && t[l] > a / 2 && (f = ue, d = h);
      }
    }
    n[a][l] = {
      str: s[c],
      val: d,
      type: f
    }, i < n[a][l].val.length + 3 && (i = n[a][l].val.length + 3);
  }
  return n[a][++l] = {
    type: ge
  }, [l + 1, i];
}
var Jn = { encode: Vn, decode: Qn }, jn = nn, ea = Ot, ta = Vt, ra = Tn, na = Jn;
function aa(n, e) {
  jn.decode(n).copy(e, 0, 0);
}
function ia(n, e) {
  ea.decode(n).copy(e, 0, 0);
}
function oa(n, e) {
  ta.decode(n).copy(e, 0, 0);
}
function sa(n, e) {
  var t = new Array();
  ra.decode(n, t).copy(e, 0, 0);
}
function fa(n, e) {
  var t = na.decode(n, 0, "\0");
  _.Buffer.from(t, "binary").copy(e, 0, 0);
}
var ca = {
  r4x8_uncompress: aa,
  r4x16_uncompress: ia,
  arith_uncompress: oa,
  fqzcomp_uncompress: sa,
  tok3_uncompress: fa
};
const Ae = /* @__PURE__ */ Bt(ca);
var tr = { exports: {} }, rr = { exports: {} };
(function() {
  var n = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", e = {
    // Bit-wise rotation left
    rotl: function(t, r) {
      return t << r | t >>> 32 - r;
    },
    // Bit-wise rotation right
    rotr: function(t, r) {
      return t << 32 - r | t >>> r;
    },
    // Swap big-endian to little-endian and vice versa
    endian: function(t) {
      if (t.constructor == Number)
        return e.rotl(t, 8) & 16711935 | e.rotl(t, 24) & 4278255360;
      for (var r = 0; r < t.length; r++)
        t[r] = e.endian(t[r]);
      return t;
    },
    // Generate an array of any length of random bytes
    randomBytes: function(t) {
      for (var r = []; t > 0; t--)
        r.push(Math.floor(Math.random() * 256));
      return r;
    },
    // Convert a byte array to big-endian 32-bit words
    bytesToWords: function(t) {
      for (var r = [], a = 0, i = 0; a < t.length; a++, i += 8)
        r[i >>> 5] |= t[a] << 24 - i % 32;
      return r;
    },
    // Convert big-endian 32-bit words to a byte array
    wordsToBytes: function(t) {
      for (var r = [], a = 0; a < t.length * 32; a += 8)
        r.push(t[a >>> 5] >>> 24 - a % 32 & 255);
      return r;
    },
    // Convert a byte array to a hex string
    bytesToHex: function(t) {
      for (var r = [], a = 0; a < t.length; a++)
        r.push((t[a] >>> 4).toString(16)), r.push((t[a] & 15).toString(16));
      return r.join("");
    },
    // Convert a hex string to a byte array
    hexToBytes: function(t) {
      for (var r = [], a = 0; a < t.length; a += 2)
        r.push(parseInt(t.substr(a, 2), 16));
      return r;
    },
    // Convert a byte array to a base-64 string
    bytesToBase64: function(t) {
      for (var r = [], a = 0; a < t.length; a += 3)
        for (var i = t[a] << 16 | t[a + 1] << 8 | t[a + 2], o = 0; o < 4; o++)
          a * 8 + o * 6 <= t.length * 8 ? r.push(n.charAt(i >>> 6 * (3 - o) & 63)) : r.push("=");
      return r.join("");
    },
    // Convert a base-64 string to a byte array
    base64ToBytes: function(t) {
      t = t.replace(/[^A-Z0-9+\/]/ig, "");
      for (var r = [], a = 0, i = 0; a < t.length; i = ++a % 4)
        i != 0 && r.push((n.indexOf(t.charAt(a - 1)) & Math.pow(2, -2 * i + 8) - 1) << i * 2 | n.indexOf(t.charAt(a)) >>> 6 - i * 2);
      return r;
    }
  };
  rr.exports = e;
})();
var la = rr.exports, et = {
  // UTF-8 encoding
  utf8: {
    // Convert a string to a byte array
    stringToBytes: function(n) {
      return et.bin.stringToBytes(unescape(encodeURIComponent(n)));
    },
    // Convert a byte array to a string
    bytesToString: function(n) {
      return decodeURIComponent(escape(et.bin.bytesToString(n)));
    }
  },
  // Binary encoding
  bin: {
    // Convert a string to a byte array
    stringToBytes: function(n) {
      for (var e = [], t = 0; t < n.length; t++)
        e.push(n.charCodeAt(t) & 255);
      return e;
    },
    // Convert a byte array to a string
    bytesToString: function(n) {
      for (var e = [], t = 0; t < n.length; t++)
        e.push(String.fromCharCode(n[t]));
      return e.join("");
    }
  }
}, St = et;
/*!
 * Determine if an object is a Buffer
 *
 * @author   Feross Aboukhadijeh <https://feross.org>
 * @license  MIT
 */
var da = function(n) {
  return n != null && (nr(n) || ha(n) || !!n._isBuffer);
};
function nr(n) {
  return !!n.constructor && typeof n.constructor.isBuffer == "function" && n.constructor.isBuffer(n);
}
function ha(n) {
  return typeof n.readFloatLE == "function" && typeof n.slice == "function" && nr(n.slice(0, 0));
}
(function() {
  var n = la, e = St.utf8, t = da, r = St.bin, a = function(i, o) {
    i.constructor == String ? o && o.encoding === "binary" ? i = r.stringToBytes(i) : i = e.stringToBytes(i) : t(i) ? i = Array.prototype.slice.call(i, 0) : !Array.isArray(i) && i.constructor !== Uint8Array && (i = i.toString());
    for (var s = n.bytesToWords(i), c = i.length * 8, l = 1732584193, f = -271733879, d = -1732584194, h = 271733878, u = 0; u < s.length; u++)
      s[u] = (s[u] << 8 | s[u] >>> 24) & 16711935 | (s[u] << 24 | s[u] >>> 8) & 4278255360;
    s[c >>> 5] |= 128 << c % 32, s[(c + 64 >>> 9 << 4) + 14] = c;
    for (var g = a._ff, m = a._gg, v = a._hh, w = a._ii, u = 0; u < s.length; u += 16) {
      var A = l, M = f, y = d, p = h;
      l = g(l, f, d, h, s[u + 0], 7, -680876936), h = g(h, l, f, d, s[u + 1], 12, -389564586), d = g(d, h, l, f, s[u + 2], 17, 606105819), f = g(f, d, h, l, s[u + 3], 22, -1044525330), l = g(l, f, d, h, s[u + 4], 7, -176418897), h = g(h, l, f, d, s[u + 5], 12, 1200080426), d = g(d, h, l, f, s[u + 6], 17, -1473231341), f = g(f, d, h, l, s[u + 7], 22, -45705983), l = g(l, f, d, h, s[u + 8], 7, 1770035416), h = g(h, l, f, d, s[u + 9], 12, -1958414417), d = g(d, h, l, f, s[u + 10], 17, -42063), f = g(f, d, h, l, s[u + 11], 22, -1990404162), l = g(l, f, d, h, s[u + 12], 7, 1804603682), h = g(h, l, f, d, s[u + 13], 12, -40341101), d = g(d, h, l, f, s[u + 14], 17, -1502002290), f = g(f, d, h, l, s[u + 15], 22, 1236535329), l = m(l, f, d, h, s[u + 1], 5, -165796510), h = m(h, l, f, d, s[u + 6], 9, -1069501632), d = m(d, h, l, f, s[u + 11], 14, 643717713), f = m(f, d, h, l, s[u + 0], 20, -373897302), l = m(l, f, d, h, s[u + 5], 5, -701558691), h = m(h, l, f, d, s[u + 10], 9, 38016083), d = m(d, h, l, f, s[u + 15], 14, -660478335), f = m(f, d, h, l, s[u + 4], 20, -405537848), l = m(l, f, d, h, s[u + 9], 5, 568446438), h = m(h, l, f, d, s[u + 14], 9, -1019803690), d = m(d, h, l, f, s[u + 3], 14, -187363961), f = m(f, d, h, l, s[u + 8], 20, 1163531501), l = m(l, f, d, h, s[u + 13], 5, -1444681467), h = m(h, l, f, d, s[u + 2], 9, -51403784), d = m(d, h, l, f, s[u + 7], 14, 1735328473), f = m(f, d, h, l, s[u + 12], 20, -1926607734), l = v(l, f, d, h, s[u + 5], 4, -378558), h = v(h, l, f, d, s[u + 8], 11, -2022574463), d = v(d, h, l, f, s[u + 11], 16, 1839030562), f = v(f, d, h, l, s[u + 14], 23, -35309556), l = v(l, f, d, h, s[u + 1], 4, -1530992060), h = v(h, l, f, d, s[u + 4], 11, 1272893353), d = v(d, h, l, f, s[u + 7], 16, -155497632), f = v(f, d, h, l, s[u + 10], 23, -1094730640), l = v(l, f, d, h, s[u + 13], 4, 681279174), h = v(h, l, f, d, s[u + 0], 11, -358537222), d = v(d, h, l, f, s[u + 3], 16, -722521979), f = v(f, d, h, l, s[u + 6], 23, 76029189), l = v(l, f, d, h, s[u + 9], 4, -640364487), h = v(h, l, f, d, s[u + 12], 11, -421815835), d = v(d, h, l, f, s[u + 15], 16, 530742520), f = v(f, d, h, l, s[u + 2], 23, -995338651), l = w(l, f, d, h, s[u + 0], 6, -198630844), h = w(h, l, f, d, s[u + 7], 10, 1126891415), d = w(d, h, l, f, s[u + 14], 15, -1416354905), f = w(f, d, h, l, s[u + 5], 21, -57434055), l = w(l, f, d, h, s[u + 12], 6, 1700485571), h = w(h, l, f, d, s[u + 3], 10, -1894986606), d = w(d, h, l, f, s[u + 10], 15, -1051523), f = w(f, d, h, l, s[u + 1], 21, -2054922799), l = w(l, f, d, h, s[u + 8], 6, 1873313359), h = w(h, l, f, d, s[u + 15], 10, -30611744), d = w(d, h, l, f, s[u + 6], 15, -1560198380), f = w(f, d, h, l, s[u + 13], 21, 1309151649), l = w(l, f, d, h, s[u + 4], 6, -145523070), h = w(h, l, f, d, s[u + 11], 10, -1120210379), d = w(d, h, l, f, s[u + 2], 15, 718787259), f = w(f, d, h, l, s[u + 9], 21, -343485551), l = l + A >>> 0, f = f + M >>> 0, d = d + y >>> 0, h = h + p >>> 0;
    }
    return n.endian([l, f, d, h]);
  };
  a._ff = function(i, o, s, c, l, f, d) {
    var h = i + (o & s | ~o & c) + (l >>> 0) + d;
    return (h << f | h >>> 32 - f) + o;
  }, a._gg = function(i, o, s, c, l, f, d) {
    var h = i + (o & c | s & ~c) + (l >>> 0) + d;
    return (h << f | h >>> 32 - f) + o;
  }, a._hh = function(i, o, s, c, l, f, d) {
    var h = i + (o ^ s ^ c) + (l >>> 0) + d;
    return (h << f | h >>> 32 - f) + o;
  }, a._ii = function(i, o, s, c, l, f, d) {
    var h = i + (s ^ (o | ~c)) + (l >>> 0) + d;
    return (h << f | h >>> 32 - f) + o;
  }, a._blocksize = 16, a._digestsize = 16, tr.exports = function(i, o) {
    if (i == null)
      throw new Error("Illegal argument " + i);
    var s = n.wordsToBytes(a(i, o));
    return o && o.asBytes ? s : o && o.asString ? r.bytesToString(s) : n.bytesToHex(s);
  };
})();
var ua = tr.exports;
const ga = /* @__PURE__ */ Bt(ua);
class me extends Error {
}
function X(n, e, t) {
  let r = 0;
  if (e.bytePosition + (7 - e.bitPosition + t) / 8 > n.length)
    throw new me("read error during decoding. the file seems to be truncated.");
  for (let a = t; a; a--)
    r <<= 1, r |= n[e.bytePosition] >> e.bitPosition & 1, e.bitPosition -= 1, e.bitPosition < 0 && (e.bytePosition += 1), e.bitPosition &= 7;
  return r;
}
function va(n) {
  return n & -128 ? n & -16384 ? n & -2097152 ? n & -268435456 ? 5 : 4 : 3 : 2 : 1;
}
function ma(n, e) {
  let t = e;
  const r = n[t];
  let a;
  if (r < 128 ? (a = r, t = t + 1) : r < 192 ? (a = (r << 8 | n[t + 1]) & 16383, t = t + 2) : r < 224 ? (a = (r << 16 | n[t + 1] << 8 | n[t + 2]) & 2097151, t = t + 3) : r < 240 ? (a = (r << 24 | n[t + 1] << 16 | n[t + 2] << 8 | n[t + 3]) & 268435455, t = t + 4) : (a = (r & 15) << 28 | n[t + 1] << 20 | n[t + 2] << 12 | n[t + 3] << 4 | n[t + 4] & 15, t = t + 5), t > n.length)
    throw new me("Attempted to read beyond end of buffer; this file seems truncated.");
  return [a, t - e];
}
function Y(n, e, t = 0, r = 0) {
  const { offset: a, result: i } = e.parse(n);
  return {
    ...i,
    _endPosition: a + r,
    _size: a - t
  };
}
function Ie(n, e) {
  const t = n.prototype[e], r = `_memo_${e}`;
  n.prototype[e] = function() {
    if (!(r in this)) {
      const i = t.call(this);
      this[r] = i, Promise.resolve(i).catch(() => {
        delete this[r];
      });
    }
    return this[r];
  };
}
function wa(n) {
  return ga(n.toUpperCase().replace(/[^\x21-\x7e]/g, ""));
}
const ya = {
  CRAM_FLAG_PRESERVE_QUAL_SCORES: 1,
  CRAM_FLAG_DETACHED: 2,
  CRAM_FLAG_MATE_DOWNSTREAM: 4,
  CRAM_FLAG_NO_SEQ: 8,
  CRAM_FLAG_MASK: 16 - 1,
  // mate read is reversed
  CRAM_M_REVERSE: 1,
  // mated read is unmapped
  CRAM_M_UNMAP: 2,
  //  the read is paired in sequencing, no matter whether it is mapped in a pair
  BAM_FPAIRED: 1,
  //  the read is mapped in a proper pair
  BAM_FPROPER_PAIR: 2,
  //  the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
  BAM_FUNMAP: 4,
  //  the mate is unmapped
  BAM_FMUNMAP: 8,
  //  the read is mapped to the reverse strand
  BAM_FREVERSE: 16,
  //  the mate is mapped to the reverse strand
  BAM_FMREVERSE: 32,
  //  this is read1
  BAM_FREAD1: 64,
  //  this is read2
  BAM_FREAD2: 128,
  //  not primary alignment
  BAM_FSECONDARY: 256,
  //  QC failure
  BAM_FQCFAIL: 512,
  //  optical or PCR duplicate
  BAM_FDUP: 1024,
  //  supplementary alignment
  BAM_FSUPPLEMENTARY: 2048,
  BAM_CMATCH: 0,
  BAM_CINS: 1,
  BAM_CDEL: 2,
  BAM_CREF_SKIP: 3,
  BAM_CSOFT_CLIP: 4,
  BAM_CHARD_CLIP: 5,
  BAM_CPAD: 6,
  BAM_CEQUAL: 7,
  BAM_CDIFF: 8,
  BAM_CBACK: 9,
  BAM_CIGAR_STR: "MIDNSHP:XB",
  BAM_CIGAR_SHIFT: 4,
  BAM_CIGAR_MASK: 15,
  BAM_CIGAR_TYPE: 246183
}, P = ya;
function pa(n, e) {
  if (!n.lengthOnRef && !n.readLength || n.isUnknownBases())
    return null;
  const t = n.alignmentStart - e.start;
  if (!n.readFeatures)
    return e.seq.substr(t, n.lengthOnRef).toUpperCase();
  let r = "", a = t, i = 0;
  for (; r.length < n.readLength; )
    if (i < n.readFeatures.length) {
      const o = n.readFeatures[i];
      if (o.code === "Q" || o.code === "q")
        i += 1;
      else if (o.pos === r.length + 1)
        if (i += 1, o.code === "b") {
          const s = o.data;
          r += s, a += s.length;
        } else
          o.code === "B" ? (r += o.data[0], a += 1) : o.code === "X" ? (r += o.sub, a += 1) : o.code === "I" ? r += o.data : o.code === "D" ? a += o.data : o.code === "i" ? r += o.data : o.code === "N" ? a += o.data : o.code === "S" ? r += o.data : o.code === "P" || o.code;
      else if (i < n.readFeatures.length) {
        const s = e.seq.substr(a, n.readFeatures[i].pos - r.length - 1);
        r += s, a += s.length;
      }
    } else {
      const o = e.seq.substr(a, n.readLength - r.length);
      r += o, a += o.length;
    }
  return r.toUpperCase();
}
const _a = {
  a: 0,
  A: 0,
  c: 1,
  C: 1,
  g: 2,
  G: 2,
  t: 3,
  T: 3,
  n: 4,
  N: 4
};
function Sa(n, e, t, r) {
  if (!e)
    return;
  const a = r.refPos - e.start, i = e.seq.charAt(a);
  i && (r.ref = i);
  let o = _a[i];
  o === void 0 && (o = 4);
  const c = t.substitutionMatrix[o][r.data];
  c && (r.sub = c);
}
const Aa = [
  [1, "Paired"],
  [2, "ProperlyPaired"],
  [4, "SegmentUnmapped"],
  [8, "MateUnmapped"],
  [16, "ReverseComplemented"],
  //  the mate is mapped to the reverse strand
  [32, "MateReverseComplemented"],
  //  this is read1
  [64, "Read1"],
  //  this is read2
  [128, "Read2"],
  //  not primary alignment
  [256, "Secondary"],
  //  QC failure
  [512, "FailedQc"],
  //  optical or PCR duplicate
  [1024, "Duplicate"],
  //  supplementary alignment
  [2048, "Supplementary"]
], ba = [
  [1, "PreservingQualityScores"],
  [2, "Detached"],
  [4, "WithMateDownstream"],
  [8, "DecodeSequenceAsStar"]
], Ca = [
  [1, "OnNegativeStrand"],
  [2, "Unmapped"]
];
function ct(n) {
  const e = {};
  for (const [t, r] of n)
    e["is" + r] = (a) => !!(a & t), e["set" + r] = (a) => a | t;
  return e;
}
const De = ct(Aa), ae = ct(ba), At = ct(Ca);
class Ba {
  constructor({ flags: e, cramFlags: t, readLength: r, mappingQuality: a, lengthOnRef: i, qualityScores: o, mateRecordNumber: s, readBases: c, readFeatures: l, mateToUse: f, readGroupId: d, readName: h, sequenceId: u, uniqueId: g, templateSize: m, alignmentStart: v, tags: w }) {
    this.flags = e, this.cramFlags = t, this.readLength = r, this.mappingQuality = a, this.lengthOnRef = i, this.qualityScores = o, c && (this.readBases = c), this.readGroupId = d, this.readName = h, this.sequenceId = u, this.uniqueId = g, this.templateSize = m, this.alignmentStart = v, this.tags = w, l && (this.readFeatures = l), f && (this.mate = {
      flags: f.mateFlags,
      readName: f.mateReadName,
      sequenceId: f.mateSequenceId,
      alignmentStart: f.mateAlignmentStart
    }), s && (this.mateRecordNumber = s);
  }
  /**
   * @returns {boolean} true if the read is paired, regardless of whether both segments are mapped
   */
  isPaired() {
    return !!(this.flags & P.BAM_FPAIRED);
  }
  /** @returns {boolean} true if the read is paired, and both segments are mapped */
  isProperlyPaired() {
    return !!(this.flags & P.BAM_FPROPER_PAIR);
  }
  /** @returns {boolean} true if the read itself is unmapped; conflictive with isProperlyPaired */
  isSegmentUnmapped() {
    return !!(this.flags & P.BAM_FUNMAP);
  }
  /** @returns {boolean} true if the read itself is unmapped; conflictive with isProperlyPaired */
  isMateUnmapped() {
    return !!(this.flags & P.BAM_FMUNMAP);
  }
  /** @returns {boolean} true if the read is mapped to the reverse strand */
  isReverseComplemented() {
    return !!(this.flags & P.BAM_FREVERSE);
  }
  /** @returns {boolean} true if the mate is mapped to the reverse strand */
  isMateReverseComplemented() {
    return !!(this.flags & P.BAM_FMREVERSE);
  }
  /** @returns {boolean} true if this is read number 1 in a pair */
  isRead1() {
    return !!(this.flags & P.BAM_FREAD1);
  }
  /** @returns {boolean} true if this is read number 2 in a pair */
  isRead2() {
    return !!(this.flags & P.BAM_FREAD2);
  }
  /** @returns {boolean} true if this is a secondary alignment */
  isSecondary() {
    return !!(this.flags & P.BAM_FSECONDARY);
  }
  /** @returns {boolean} true if this read has failed QC checks */
  isFailedQc() {
    return !!(this.flags & P.BAM_FQCFAIL);
  }
  /** @returns {boolean} true if the read is an optical or PCR duplicate */
  isDuplicate() {
    return !!(this.flags & P.BAM_FDUP);
  }
  /** @returns {boolean} true if this is a supplementary alignment */
  isSupplementary() {
    return !!(this.flags & P.BAM_FSUPPLEMENTARY);
  }
  /**
   * @returns {boolean} true if the read is detached
   */
  isDetached() {
    return !!(this.cramFlags & P.CRAM_FLAG_DETACHED);
  }
  /** @returns {boolean} true if the read has a mate in this same CRAM segment */
  hasMateDownStream() {
    return !!(this.cramFlags & P.CRAM_FLAG_MATE_DOWNSTREAM);
  }
  /** @returns {boolean} true if the read contains qual scores */
  isPreservingQualityScores() {
    return !!(this.cramFlags & P.CRAM_FLAG_PRESERVE_QUAL_SCORES);
  }
  /** @returns {boolean} true if the read has no sequence bases */
  isUnknownBases() {
    return !!(this.cramFlags & P.CRAM_FLAG_NO_SEQ);
  }
  /**
   * Get the original sequence of this read.
   * @returns {String} sequence basepairs
   */
  getReadBases() {
    if (!this.readBases && this._refRegion) {
      const e = pa(this, this._refRegion);
      e && (this.readBases = e);
    }
    return this.readBases;
  }
  /**
   * Get the pair orientation of a paired read. Adapted from igv.js
   * @returns {String} of paired orientatin
   */
  getPairOrientation() {
    if (!this.isSegmentUnmapped() && this.isPaired() && !this.isMateUnmapped() && this.mate && this.sequenceId === this.mate.sequenceId) {
      const e = this.isReverseComplemented() ? "R" : "F", t = this.isMateReverseComplemented() ? "R" : "F";
      let r = " ", a = " ";
      this.isRead1() ? (r = "1", a = "2") : this.isRead2() && (r = "2", a = "1");
      const i = [];
      let o = this.templateLength || this.templateSize;
      if (o === void 0)
        throw new Error("One of templateSize and templateLength must be set");
      return this.alignmentStart > this.mate.alignmentStart && o > 0 && (o = -o), o > 0 ? (i[0] = e, i[1] = r, i[2] = t, i[3] = a) : (i[2] = e, i[3] = r, i[0] = t, i[1] = a), i.join("");
    }
    return null;
  }
  /**
   * Annotates this feature with the given reference sequence basepair
   * information. This will add a `sub` and a `ref` item to base
   * subsitution read features given the actual substituted and reference
   * base pairs, and will make the `getReadSequence()` method work.
   *
   * @param {object} refRegion
   * @param {number} refRegion.start
   * @param {number} refRegion.end
   * @param {string} refRegion.seq
   * @param {CramContainerCompressionScheme} compressionScheme
   * @returns {undefined} nothing
   */
  addReferenceSequence(e, t) {
    this.readFeatures && this.readFeatures.forEach((r) => {
      r.code === "X" && Sa(this, e, t, r);
    }), !this.readBases && e.start <= this.alignmentStart && e.end >= this.alignmentStart + (this.lengthOnRef || this.readLength) - 1 && (this._refRegion = e);
  }
  toJSON() {
    const e = {};
    return Object.keys(this).forEach((t) => {
      t.charAt(0) !== "_" && (e[t] = this[t]);
    }), e.readBases = this.getReadBases(), e;
  }
}
function xe(n) {
  let e = "";
  for (let t = 0; t < n.length && n[t] !== 0; t++)
    e += String.fromCharCode(n[t]);
  return e;
}
function Ra(n) {
  const e = String.fromCharCode(n[0]), t = Int32Array.from(n.slice(1))[0], r = new Array(t);
  if (n = n.slice(5), e === "c") {
    const a = new Int8Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else if (e === "C") {
    const a = new Uint8Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else if (e === "s") {
    const a = new Int16Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else if (e === "S") {
    const a = new Uint16Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else if (e === "i") {
    const a = new Int32Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else if (e === "I") {
    const a = new Uint32Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else if (e === "f") {
    const a = new Float32Array(n.buffer);
    for (let i = 0; i < t; i += 1)
      r[i] = a[i];
  } else
    throw new Error("unknown type: " + e);
  return r;
}
function Ea(n, e) {
  if (n === "Z")
    return xe(e);
  if (n === "A")
    return String.fromCharCode(e[0]);
  if (n === "I")
    return yr.fromBytesLE(e).toNumber();
  if (n === "i")
    return new Int32Array(e.buffer)[0];
  if (n === "s")
    return new Int16Array(e.buffer)[0];
  if (n === "S")
    return new Uint16Array(e.buffer)[0];
  if (n === "c")
    return new Int8Array(e.buffer)[0];
  if (n === "C")
    return e[0];
  if (n === "f")
    return new Float32Array(e.buffer)[0];
  if (n === "H")
    return Number.parseInt(xe(e).replace(/^0x/, ""), 16);
  if (n === "B")
    return Ra(e);
  throw new C(`Unrecognized tag type ${n}`);
}
function xa(n, e, t, r, a) {
  let i = 0, o = n - 1;
  const s = new Array(e);
  function c([l, f]) {
    const d = t(f);
    if (l === "character")
      return String.fromCharCode(d);
    if (l === "string") {
      let h = "";
      for (let u = 0; u < d.byteLength; u++)
        h += String.fromCharCode(d[u]);
      return h;
    }
    return l === "numArray" ? d.toArray() : d;
  }
  for (let l = 0; l < e; l += 1) {
    const f = String.fromCharCode(t("FC")), d = t("FP"), h = {
      B: ["character", "BA"],
      S: ["string", a > 1 ? "SC" : "IN"],
      X: ["number", "BS"],
      D: ["number", "DL"],
      I: ["string", "IN"],
      i: ["character", "BA"],
      b: ["string", "BB"],
      q: ["numArray", "QQ"],
      Q: ["number", "QS"],
      H: ["number", "HC"],
      P: ["number", "PD"],
      N: ["number", "RS"]
    }[f];
    if (!h)
      throw new C(`invalid read feature code "${f}"`);
    let u = c(h);
    const g = { B: ["number", "QS"] }[f];
    g && (u = [u, c(g)]), i += d;
    const m = i;
    o += d;
    const v = o;
    f === "D" || f === "N" ? o += u : f === "I" || f === "S" ? o -= u.length : f === "i" && (o -= 1), s[l] = { code: f, pos: m, refPos: v, data: u };
  }
  return s;
}
function Ia(n, e, t, r, a, i, o, s, c) {
  let l = e("BF");
  const f = e("CF");
  if (!ie(r.parsedContent))
    throw new Error();
  let d;
  s > 1 && r.parsedContent.refSeqId === -2 ? d = e("RI") : d = r.parsedContent.refSeqId;
  const h = e("RL");
  let u = e("AP");
  t.APdelta && (u = u + o.lastAlignmentStart), o.lastAlignmentStart = u;
  const g = e("RG");
  let m;
  t.readNamesIncluded && (m = xe(e("RN")));
  let v, w, A;
  if (ae.isDetached(f)) {
    const F = e("MF");
    let b;
    t.readNamesIncluded || (b = xe(e("RN")), m = b);
    const D = e("NS"), $ = e("NP");
    (F || D > -1) && (v = {
      mateFlags: F,
      mateSequenceId: D,
      mateAlignmentStart: $,
      mateReadName: b
    }), w = e("TS"), At.isUnmapped(F) && (l = De.setMateUnmapped(l)), At.isOnNegativeStrand(F) && (l = De.setMateReverseComplemented(l));
  } else
    ae.isWithMateDownstream(f) && (A = e("NF") + c + 1);
  const M = e("TL");
  if (M < 0)
    throw new C("invalid TL index");
  const y = {}, p = t.getTagNames(M), B = p.length;
  for (let F = 0; F < B; F += 1) {
    const b = p[F], D = b.substr(0, 2), $ = b.substr(2, 1), N = t.getCodecForTag(b);
    if (!N)
      throw new C(`no codec defined for auxiliary tag ${b}`);
    const T = N.decode(n, a, i, o);
    y[D] = Ea($, T);
  }
  let k, I, L, R, W;
  if (De.isSegmentUnmapped(l))
    if (ae.isDecodeSequenceAsStar(f))
      W = null, R = null;
    else {
      const F = new Array(h);
      for (let b = 0; b < F.length; b += 1)
        F[b] = e("BA");
      if (W = String.fromCharCode(...F), ae.isPreservingQualityScores(f)) {
        R = new Array(h);
        for (let b = 0; b < F.length; b += 1)
          R[b] = e("QS");
      }
    }
  else {
    const F = e("FN");
    if (F && (k = xa(u, F, e, t, s)), I = h, k)
      for (const { code: b, data: D } of k)
        b === "D" || b === "N" ? I += D : b === "I" || b === "S" ? I = I - D.length : b === "i" && (I = I - 1);
    if (Number.isNaN(I) && (console.warn(`${m || `${d}:${u}`} record has invalid read features`), I = h), L = e("MQ"), ae.isPreservingQualityScores(f)) {
      R = new Array(h);
      for (let b = 0; b < R.length; b++)
        R[b] = e("QS");
    }
  }
  return {
    readLength: h,
    sequenceId: d,
    cramFlags: f,
    flags: l,
    alignmentStart: u,
    readGroupId: g,
    readName: m,
    mateToUse: v,
    templateSize: w,
    mateRecordNumber: A,
    readFeatures: k,
    lengthOnRef: I,
    mappingQuality: L,
    qualityScores: R,
    readBases: W,
    tags: y
  };
}
function Ma(n, e, t) {
  function r(c) {
    const l = [c];
    if (c.mateRecordNumber !== void 0 && c.mateRecordNumber >= 0) {
      const f = n[c.mateRecordNumber];
      if (!f)
        throw new C("intra-slice mate record not found, this file seems malformed");
      l.push(...r(f));
    }
    return l;
  }
  const a = r(t), i = a.map((c) => c.alignmentStart), o = a.map((c) => c.alignmentStart + c.readLength - 1), s = Math.max(...o) - Math.min(...i) + 1;
  s >= 0 && a.forEach((c) => {
    if (c.templateLength !== void 0)
      throw new C("mate pair group has some members that have template lengths already, this file seems malformed");
    c.templateLength = s;
  });
}
function Fa(n, e) {
  const t = Math.min(n.alignmentStart, e.alignmentStart), a = Math.max(n.alignmentStart + n.readLength - 1, e.alignmentStart + e.readLength - 1) - t + 1;
  n.templateLength = a, e.templateLength = a;
}
function qa(n, e, t, r) {
  if (!r)
    throw new C("could not resolve intra-slice mate pairs, file seems truncated or malformed");
  const a = !!(r.mate || r.mateRecordNumber !== void 0 && r.mateRecordNumber !== e);
  t.readName || (t.readName = String(t.uniqueId), r.readName = t.readName), t.mate = {
    sequenceId: r.sequenceId,
    alignmentStart: r.alignmentStart,
    uniqueId: r.uniqueId
  }, r.readName && (t.mate.readName = r.readName), !r.mate && r.mateRecordNumber === void 0 && (r.mate = {
    sequenceId: t.sequenceId,
    alignmentStart: t.alignmentStart,
    uniqueId: t.uniqueId
  }, t.readName && (r.mate.readName = t.readName)), t.flags |= P.BAM_FPAIRED, r.flags & P.BAM_FUNMAP && (t.flags |= P.BAM_FMUNMAP), t.flags & P.BAM_FUNMAP && (r.flags |= P.BAM_FMUNMAP), r.flags & P.BAM_FREVERSE && (t.flags |= P.BAM_FMREVERSE), t.flags & P.BAM_FREVERSE && (r.flags |= P.BAM_FMREVERSE), t.templateLength === void 0 && (a ? Ma(n, e, t) : Fa(t, r)), delete t.mateRecordNumber;
}
class ar {
  constructor(e, t, r) {
    this.container = e, this.containerPosition = t, this.file = e.file;
  }
  // memoize
  async getHeader() {
    const e = await this.file.getSectionParsers(), t = await this.container.getHeader(), r = await this.file.readBlock(t._endPosition + this.containerPosition);
    if (r === void 0)
      throw new Error();
    if (r.contentType === "MAPPED_SLICE_HEADER") {
      const a = Y(r.content, e.cramMappedSliceHeader.parser, 0, t._endPosition);
      return { ...r, parsedContent: a };
    } else if (r.contentType === "UNMAPPED_SLICE_HEADER") {
      const a = Y(r.content, e.cramUnmappedSliceHeader.parser, 0, t._endPosition);
      return { ...r, parsedContent: a };
    } else
      throw new C(`error reading slice header block, invalid content type ${r.contentType}`);
  }
  // memoize
  async getBlocks() {
    const e = await this.getHeader();
    let t = e._endPosition;
    const r = new Array(e.parsedContent.numBlocks);
    for (let a = 0; a < r.length; a += 1) {
      const i = await this.file.readBlock(t);
      if (i === void 0)
        throw new Error();
      r[a] = i, t = r[a]._endPosition;
    }
    return r;
  }
  // no memoize
  async getCoreDataBlock() {
    return (await this.getBlocks())[0];
  }
  // memoize
  async _getBlocksContentIdIndex() {
    const e = await this.getBlocks(), t = {};
    return e.forEach((r) => {
      r.contentType === "EXTERNAL_DATA" && (t[r.contentId] = r);
    }), t;
  }
  async getBlockByContentId(e) {
    return (await this._getBlocksContentIdIndex())[e];
  }
  async getReferenceRegion() {
    const e = (await this.getHeader()).parsedContent;
    if (!ie(e))
      throw new Error();
    if (e.refSeqId < 0)
      return;
    const t = await this.container.getCompressionScheme();
    if (t === void 0)
      throw new Error();
    if (e.refBaseBlockId >= 0) {
      const r = await this.getBlockByContentId(e.refBaseBlockId);
      if (!r)
        throw new C("embedded reference specified, but reference block does not exist");
      return {
        seq: r.data.toString("utf8"),
        start: e.refSeqStart,
        end: e.refSeqStart + e.refSeqSpan - 1,
        span: e.refSeqSpan
      };
    }
    if (t.referenceRequired || this.file.fetchReferenceSequenceCallback) {
      if (!this.file.fetchReferenceSequenceCallback)
        throw new Error("reference sequence not embedded, and seqFetch callback not provided, cannot fetch reference sequence");
      const r = await this.file.fetchReferenceSequenceCallback(e.refSeqId, e.refSeqStart, e.refSeqStart + e.refSeqSpan - 1);
      if (r.length !== e.refSeqSpan)
        throw new Sr("seqFetch callback returned a reference sequence of the wrong length");
      return {
        seq: r,
        start: e.refSeqStart,
        end: e.refSeqStart + e.refSeqSpan - 1,
        span: e.refSeqSpan
      };
    }
  }
  getAllRecords() {
    return this.getRecords(() => !0);
  }
  async _fetchRecords() {
    const { majorVersion: e } = await this.file.getDefinition(), t = await this.container.getCompressionScheme();
    if (t === void 0)
      throw new Error();
    const r = await this.getHeader();
    if (r === void 0)
      throw new Error();
    const a = await this._getBlocksContentIdIndex();
    if (e > 1 && this.file.options.checkSequenceMD5 && ie(r.parsedContent) && r.parsedContent.refSeqId >= 0 && r.parsedContent.md5.join("") !== "0000000000000000") {
      const l = await this.getReferenceRegion();
      if (l) {
        const { seq: f, start: d, end: h } = l, u = wa(f), g = r.parsedContent.md5.map((m) => (m < 16 ? "0" : "") + m.toString(16)).join("");
        if (u !== g)
          throw new C(`MD5 checksum reference mismatch for ref ${r.parsedContent.refSeqId} pos ${d}..${h}. recorded MD5: ${g}, calculated MD5: ${u}`);
      }
    }
    const i = await this.getCoreDataBlock(), o = {
      lastAlignmentStart: ie(r.parsedContent) ? r.parsedContent.refSeqStart : 0,
      coreBlock: { bitPosition: 7, bytePosition: 0 },
      externalBlocks: {
        map: /* @__PURE__ */ new Map(),
        getCursor(l) {
          let f = this.map.get(l);
          return f === void 0 && (f = { bitPosition: 7, bytePosition: 0 }, this.map.set(l, f)), f;
        }
      }
    }, s = (l) => {
      const f = t.getCodecForDataSeries(l);
      if (!f)
        throw new C(`no codec defined for ${l} data series`);
      return f.decode(this, i, a, o);
    };
    let c = new Array(r.parsedContent.numRecords);
    for (let l = 0; l < c.length; l += 1)
      try {
        const f = Ia(this, s, t, r, i, a, o, e, l);
        c[l] = new Ba({
          ...f,
          uniqueId: r.contentPosition + r.parsedContent.recordCounter + l + 1
        });
      } catch (f) {
        if (f instanceof me) {
          console.warn("read attempted beyond end of buffer, file seems truncated."), c = c.filter((d) => !!d);
          break;
        } else
          throw f;
      }
    for (let l = 0; l < c.length; l += 1) {
      const { mateRecordNumber: f } = c[l];
      f !== void 0 && f >= 0 && qa(c, l, c[l], c[f]);
    }
    return c;
  }
  async getRecords(e) {
    const t = this.container.filePosition + this.containerPosition;
    let r = this.file.featureCache.get(t.toString());
    r || (r = this._fetchRecords(), this.file.featureCache.set(t.toString(), r));
    const i = (await r).filter(e);
    if (i.length && this.file.fetchReferenceSequenceCallback) {
      const o = await this.getHeader();
      if (ie(o.parsedContent) && (o.parsedContent.refSeqId >= 0 || // single-ref slice
      o.parsedContent.refSeqId === -2)) {
        const s = o.parsedContent.refSeqId >= 0 ? o.parsedContent.refSeqId : void 0, c = await this.container.getCompressionScheme();
        if (c === void 0)
          throw new Error();
        const l = {};
        for (let f = 0; f < i.length; f += 1) {
          const d = s !== void 0 ? s : i[f].sequenceId;
          let h = l[d];
          h || (h = {
            id: d,
            start: i[f].alignmentStart,
            end: -1 / 0,
            seq: null
          }, l[d] = h);
          const u = i[f].alignmentStart + (i[f].lengthOnRef || i[f].readLength) - 1;
          u > h.end && (h.end = u), i[f].alignmentStart < h.start && (h.start = i[f].alignmentStart);
        }
        await Promise.all(Object.values(l).map(async (f) => {
          f.id !== -1 && f.start <= f.end && (f.seq = await this.file.fetchReferenceSequenceCallback(f.id, f.start, f.end));
        }));
        for (let f = 0; f < i.length; f += 1) {
          const d = s !== void 0 ? s : i[f].sequenceId, h = l[d];
          if (h && h.seq) {
            const u = h.seq;
            i[f].addReferenceSequence({ ...h, seq: u }, c);
          }
        }
      }
    }
    return i;
  }
}
"getHeader getBlocks _getBlocksContentIdIndex".split(" ").forEach((n) => Ie(ar, n));
class Z {
  constructor(e, t) {
    this.parameters = e, this.dataType = t;
  }
}
function Pa(n) {
  let e = n - (n >> 1) & 1431655765;
  return e = (e & 858993459) + (e >> 2 & 858993459), (e + (e >> 4) & 252645135) * 16843009 >> 24;
}
class La extends Z {
  constructor(e, t) {
    if (super(e, t), this.codes = {}, this.codeBook = {}, this.sortedByValue = [], this.sortedCodes = [], this.sortedValuesByBitCode = [], this.sortedBitCodes = [], this.sortedBitLengthsByBitCode = [], this.bitCodeToValue = [], !["byte", "int"].includes(this.dataType))
      throw new TypeError(`${this.dataType} decoding not yet implemented by HUFFMAN_INT codec`);
    this.buildCodeBook(), this.buildCodes(), this.buildCaches(), this.sortedCodes[0].bitLength === 0 && (this._decode = this._decodeZeroLengthCode);
  }
  buildCodeBook() {
    let e = new Array(this.parameters.numCodes);
    for (let t = 0; t < this.parameters.numCodes; t += 1)
      e[t] = {
        symbol: this.parameters.symbols[t],
        bitLength: this.parameters.bitLengths[t]
      };
    e = e.sort((t, r) => t.bitLength - r.bitLength || t.symbol - r.symbol), this.codeBook = {}, e.forEach((t) => {
      this.codeBook[t.bitLength] || (this.codeBook[t.bitLength] = []), this.codeBook[t.bitLength].push(t.symbol);
    });
  }
  buildCodes() {
    this.codes = {};
    let e = 0, t = -1;
    Object.entries(this.codeBook).forEach(([r, a]) => {
      const i = parseInt(r, 10);
      a.forEach((o) => {
        const s = {
          bitLength: i,
          value: o,
          bitCode: 0
        };
        t = t + 1;
        const c = i - e;
        if (t = t << c, s.bitCode = t, e = e + c, Pa(t) > i)
          throw new C("Symbol out of range");
        this.codes[o] = s;
      });
    });
  }
  buildCaches() {
    this.sortedCodes = Object.values(this.codes).sort((t, r) => t.bitLength - r.bitLength || t.bitCode - r.bitCode), this.sortedByValue = Object.values(this.codes).sort((t, r) => t.value - r.value), this.sortedValuesByBitCode = this.sortedCodes.map((t) => t.value), this.sortedBitCodes = this.sortedCodes.map((t) => t.bitCode), this.sortedBitLengthsByBitCode = this.sortedCodes.map((t) => t.bitLength);
    const e = Math.max(...this.sortedBitCodes);
    this.bitCodeToValue = new Array(e + 1).fill(-1);
    for (let t = 0; t < this.sortedBitCodes.length; t += 1)
      this.bitCodeToValue[this.sortedCodes[t].bitCode] = t;
  }
  decode(e, t, r, a) {
    return this._decode(e, t, a.coreBlock);
  }
  // _decodeNull() {
  //   return -1
  // }
  // the special case for zero-length codes
  _decodeZeroLengthCode() {
    return this.sortedCodes[0].value;
  }
  _decode(e, t, r) {
    const a = t.content;
    let i = 0, o = 0;
    for (let s = 0; s < this.sortedCodes.length; s += 1) {
      const c = this.sortedCodes[s].bitLength;
      o <<= c - i, o |= X(a, r, c - i), i = c;
      {
        const l = this.bitCodeToValue[o];
        if (l > -1 && this.sortedBitLengthsByBitCode[l] === c)
          return this.sortedValuesByBitCode[l];
        for (let f = s; this.sortedCodes[f + 1].bitLength === c && f < this.sortedCodes.length; f += 1)
          s += 1;
      }
    }
    throw new C("Huffman symbol not found.");
  }
}
class ka extends Z {
  constructor(e, t) {
    if (super(e, t), this.dataType === "int")
      this._decodeData = this._decodeInt;
    else if (this.dataType === "byte")
      this._decodeData = this._decodeByte;
    else
      throw new G(`${this.dataType} decoding not yet implemented by EXTERNAL codec`);
  }
  decode(e, t, r, a) {
    const { blockContentId: i } = this.parameters, o = r[i];
    if (!o)
      throw new C(`no block found with content ID ${i}`);
    const s = a.externalBlocks.getCursor(i);
    return this._decodeData(o, s);
  }
  _decodeInt(e, t) {
    const [r, a] = ma(e.content, t.bytePosition);
    return t.bytePosition = t.bytePosition + a, r;
  }
  _decodeByte(e, t) {
    if (t.bytePosition >= e.content.length)
      throw new me("attempted to read beyond end of block. this file seems truncated.");
    return e.content[t.bytePosition++];
  }
}
let Ua = class extends Z {
  constructor(e, t) {
    if (super(e, t), t !== "byteArray")
      throw new TypeError(`byteArrayStop codec does not support data type ${t}`);
  }
  decode(e, t, r, a) {
    const { blockContentId: i } = this.parameters, o = r[i];
    if (!o)
      throw new C(`no block found with content ID ${i}`);
    const s = a.externalBlocks.getCursor(i);
    return this._decodeByteArray(o, s);
  }
  _decodeByteArray(e, t) {
    const r = e.content, { stopByte: a } = this.parameters, i = t.bytePosition;
    let o = t.bytePosition;
    for (; r[o] !== a && o < r.length; ) {
      if (o === r.length)
        throw new me("byteArrayStop reading beyond length of data buffer?");
      o = o + 1;
    }
    return t.bytePosition = o + 1, r.subarray(i, o);
  }
};
class ir extends Z {
  constructor(e, t, r) {
    if (super(e, t), this.instantiateCodec = r, t !== "byteArray")
      throw new TypeError(`byteArrayLength does not support data type ${t}`);
  }
  decode(e, t, r, a) {
    const o = this._getLengthCodec().decode(e, t, r, a), s = this._getDataCodec(), c = new Uint8Array(o);
    for (let l = 0; l < o; l += 1)
      c[l] = s.decode(e, t, r, a);
    return c;
  }
  // memoize
  _getLengthCodec() {
    const e = this.parameters.lengthsEncoding;
    return this.instantiateCodec(e, "int");
  }
  // memoize
  _getDataCodec() {
    const e = this.parameters.valuesEncoding;
    return this.instantiateCodec(e, "byte");
  }
}
"_getLengthCodec _getDataCodec".split(" ").forEach((n) => Ie(ir, n));
class Na extends Z {
  constructor(e, t) {
    if (super(e, t), this.dataType !== "int")
      throw new G(`${this.dataType} decoding not yet implemented by BETA codec`);
  }
  decode(e, t, r, a) {
    return X(t.content, a.coreBlock, this.parameters.length) - this.parameters.offset;
  }
}
class Da extends Z {
  constructor(e, t) {
    if (super(e, t), this.dataType !== "int")
      throw new G(`${this.dataType} decoding not yet implemented by GAMMA codec`);
  }
  decode(e, t, r, a) {
    let i = 1;
    for (; X(t.content, a.coreBlock, 1) === 0; )
      i = i + 1;
    return (X(t.content, a.coreBlock, i - 1) | 1 << i - 1) - this.parameters.offset;
  }
}
class $a extends Z {
  constructor(e, t) {
    if (super(e, t), this.dataType !== "int")
      throw new G(`${this.dataType} decoding not yet implemented by SUBEXP codec`);
  }
  decode(e, t, r, a) {
    let i = 0;
    for (; X(t.content, a.coreBlock, 1); )
      i = i + 1;
    let o, s;
    if (i === 0)
      o = this.parameters.K, s = X(t.content, a.coreBlock, o);
    else {
      o = i + this.parameters.K - 1;
      const c = X(t.content, a.coreBlock, o);
      s = 1 << o | c;
    }
    return s - this.parameters.offset;
  }
}
const za = {
  1: ka,
  // 2: GolombCodec,
  3: La,
  4: ir,
  5: Ua,
  6: Na,
  7: $a,
  // 8: GolombRiceCodec,
  9: Da
};
function Wa(n) {
  return za[n];
}
function tt(n, e) {
  const t = Wa(e === "ignore" ? 0 : n.codecId);
  if (!t)
    throw new G(`no codec implemented for codec ID ${n.codecId}`);
  return new t(n.parameters, e, tt);
}
const Ta = {
  BF: "int",
  CF: "int",
  RI: "int",
  RL: "int",
  AP: "int",
  RG: "int",
  MF: "int",
  NS: "int",
  NP: "int",
  TS: "int",
  NF: "int",
  TC: "byte",
  TN: "int",
  FN: "int",
  FC: "byte",
  FP: "int",
  BS: "byte",
  IN: "byteArray",
  SC: "byteArray",
  DL: "int",
  BA: "byte",
  BB: "byteArray",
  RS: "int",
  PD: "int",
  HC: "int",
  MQ: "int",
  RN: "byteArray",
  QS: "byte",
  QQ: "byteArray",
  TL: "int"
  // TM: 'ignore',
  // TV: 'ignore',
};
function Ha(n) {
  const e = new Array(5);
  for (let t = 0; t < 5; t += 1)
    e[t] = new Array(4);
  return e[0][n[0] >> 6 & 3] = "C", e[0][n[0] >> 4 & 3] = "G", e[0][n[0] >> 2 & 3] = "T", e[0][n[0] >> 0 & 3] = "N", e[1][n[1] >> 6 & 3] = "A", e[1][n[1] >> 4 & 3] = "G", e[1][n[1] >> 2 & 3] = "T", e[1][n[1] >> 0 & 3] = "N", e[2][n[2] >> 6 & 3] = "A", e[2][n[2] >> 4 & 3] = "C", e[2][n[2] >> 2 & 3] = "T", e[2][n[2] >> 0 & 3] = "N", e[3][n[3] >> 6 & 3] = "A", e[3][n[3] >> 4 & 3] = "C", e[3][n[3] >> 2 & 3] = "G", e[3][n[3] >> 0 & 3] = "N", e[4][n[4] >> 6 & 3] = "A", e[4][n[4] >> 4 & 3] = "C", e[4][n[4] >> 2 & 3] = "G", e[4][n[4] >> 0 & 3] = "T", e;
}
class Oa {
  constructor(e) {
    this.dataSeriesCodecCache = {}, this.tagCodecCache = {}, this.tagEncoding = {}, this.readNamesIncluded = e.preservation.RN, this.APdelta = e.preservation.AP, this.referenceRequired = !!e.preservation.RR, this.tagIdsDictionary = e.preservation.TD, this.substitutionMatrix = Ha(e.preservation.SM), this.dataSeriesEncoding = e.dataSeriesEncoding, this.tagEncoding = e.tagEncoding, this.preservation = e.preservation, this._size = e._size, this._endPosition = e._endPosition;
  }
  /**
   * @param {string} tagName three-character tag name
   * @private
   */
  getCodecForTag(e) {
    if (!this.tagCodecCache[e]) {
      const t = this.tagEncoding[e];
      t && (this.tagCodecCache[e] = tt(t, "byteArray"));
    }
    return this.tagCodecCache[e];
  }
  /**
   *
   * @param {number} tagListId ID of the tag list to fetch from the tag dictionary
   * @private
   */
  getTagNames(e) {
    return this.tagIdsDictionary[e];
  }
  getCodecForDataSeries(e) {
    let t = this.dataSeriesCodecCache[e];
    if (t === void 0) {
      const r = this.dataSeriesEncoding[e];
      if (r) {
        const a = Ta[e];
        if (!a)
          throw new C(`data series name ${e} not defined in file compression header`);
        t = tt(r, a), this.dataSeriesCodecCache[e] = t;
      }
    }
    return t;
  }
  toJSON() {
    const e = {};
    return Object.keys(this).forEach((t) => {
      /Cache$/.test(t) || (e[t] = this[t]);
    }), e;
  }
}
class or {
  constructor(e, t) {
    this.file = e, this.filePosition = t;
  }
  // memoize
  getHeader() {
    return this._readContainerHeader(this.filePosition);
  }
  // memoize
  async getCompressionHeaderBlock() {
    if (!(await this.getHeader()).numRecords)
      return null;
    const t = await this.file.getSectionParsers(), r = await this.getFirstBlock();
    if (r === void 0)
      return;
    if (r.contentType !== "COMPRESSION_HEADER")
      throw new C(`invalid content type ${r.contentType} in what is supposed to be the compression header block`);
    const a = Y(r.content, t.cramCompressionHeader.parser, 0, r.contentPosition);
    return {
      ...r,
      parsedContent: a
    };
  }
  async getFirstBlock() {
    const e = await this.getHeader();
    return this.file.readBlock(e._endPosition);
  }
  // parses the compression header data into a CramContainerCompressionScheme object
  // memoize
  async getCompressionScheme() {
    const e = await this.getCompressionHeaderBlock();
    if (e)
      return new Oa(e.parsedContent);
  }
  getSlice(e, t) {
    return new ar(this, e, t);
  }
  async _readContainerHeader(e) {
    const t = await this.file.getSectionParsers(), { cramContainerHeader1: r, cramContainerHeader2: a } = t, { size: i } = await this.file.stat();
    if (e >= i)
      return;
    const o = _.Buffer.allocUnsafe(r.maxLength);
    await this.file.read(o, 0, r.maxLength, e);
    const s = Y(o, r.parser), c = va(s.numLandmarks);
    if (e + s.length >= i) {
      console.warn(`${this.file}: container header at ${e} indicates that the container has length ${s.length}, which extends beyond the length of the file. Skipping this container.`);
      return;
    }
    const l = _.Buffer.allocUnsafe(a.maxLength(s.numLandmarks));
    await this.file.read(l, 0, a.maxLength(s.numLandmarks), e + s._size - c);
    const f = Y(l, a.parser);
    return this.file.validateChecksums && f.crc32 !== void 0 && await this.file.checkCrc32(e, s._size + f._size - c - 4, f.crc32, `container header beginning at position ${e}`), Object.assign(s, f, {
      _size: s._size + f._size - c,
      _endPosition: s._size + f._size - c + e
    });
  }
}
"getHeader getCompressionHeaderBlock getCompressionScheme".split(" ").forEach((n) => Ie(or, n));
function Ga(n) {
  if (n == null)
    throw new Error("Value must not be nullish.");
  return n;
}
function Qa(n) {
  const { protocol: e, pathname: t } = pr.parse(n);
  return e === "file:" ? new Rt(unescape(Ga(t))) : new lr(n);
}
function sr(n, e, t) {
  if (t)
    return t;
  if (n)
    return Qa(n);
  if (e)
    return new Rt(e);
  throw new Error("no url, path, or filehandle provided, cannot open");
}
function bt(n) {
  const e = n.split(/\r?\n/), t = [];
  return e.forEach((r) => {
    const [a, ...i] = r.split(/\t/), o = i.map((s) => {
      const [c, l] = s.split(":", 2);
      return { tag: c, value: l };
    });
    a && t.push({ tag: a.substr(1), data: o });
  }), t;
}
function Va() {
  const n = new Uint32Array([287454020]), e = new Uint8Array(n.buffer);
  return e[0] === 68 ? 0 : e[0] === 17 ? 1 : 2;
}
class rt {
  constructor(e) {
    var t;
    if (this.file = sr(e.url, e.path, e.filehandle), this.validateChecksums = !0, this.fetchReferenceSequenceCallback = e.seqFetch, this.options = {
      checkSequenceMD5: e.checkSequenceMD5,
      cacheSize: (t = e.cacheSize) !== null && t !== void 0 ? t : 2e4
    }, this.featureCache = new Et({
      maxSize: this.options.cacheSize
    }), Va() > 0)
      throw new Error("Detected big-endian machine, may be unable to run");
  }
  // toString() {
  //   if (this.file.filename) {
  //     return this.file.filename
  //   }
  //   if (this.file.url) {
  //     return this.file.url
  //   }
  //
  //   return '(cram file)'
  // }
  // can just read this object like a filehandle
  read(e, t, r, a) {
    return this.file.read(e, t, r, a);
  }
  // can just stat this object like a filehandle
  stat() {
    return this.file.stat();
  }
  // memoized
  async getDefinition() {
    const e = _.Buffer.allocUnsafe(be.maxLength);
    await this.file.read(e, 0, be.maxLength, 0);
    const t = be.parser.parse(e).result;
    if (t.majorVersion !== 2 && t.majorVersion !== 3)
      throw new G(`CRAM version ${t.majorVersion} not supported`);
    return t;
  }
  // memoize
  async getSamHeader() {
    const e = await this.getContainerById(0);
    if (!e)
      throw new C("file contains no containers");
    const t = await e.getFirstBlock();
    if (t === void 0)
      return bt("");
    const r = t.content, a = r.readInt32LE(0), i = 4, o = r.toString("utf8", i, i + a);
    return this.header = o, bt(o);
  }
  async getHeaderText() {
    return await this.getSamHeader(), this.header;
  }
  // memoize
  async getSectionParsers() {
    const { majorVersion: e } = await this.getDefinition();
    return Hr(e);
  }
  async getContainerById(e) {
    const t = await this.getSectionParsers();
    let r = t.cramFileDefinition.maxLength;
    const { size: a } = await this.file.stat(), { cramContainerHeader1: i } = t;
    let o;
    for (let s = 0; s <= e; s += 1) {
      if (r + i.maxLength + 8 >= a)
        return;
      o = this.getContainerAtPosition(r);
      const c = await o.getHeader();
      if (!c)
        throw new C(`container ${e} not found in file`);
      if (s === 0) {
        r = c._endPosition;
        for (let l = 0; l < c.numBlocks; l += 1) {
          const f = await this.readBlock(r);
          if (f === void 0)
            return;
          r = f._endPosition;
        }
      } else
        r += c._size + c.length;
    }
    return o;
  }
  async checkCrc32(e, t, r, a) {
    const i = _.Buffer.allocUnsafe(t);
    await this.file.read(i, 0, t, e);
    const o = dr.unsigned(i);
    if (o !== r)
      throw new C(`crc mismatch in ${a}: recorded CRC32 = ${r}, but calculated CRC32 = ${o}`);
  }
  /**
   * @returns {Promise[number]} the number of containers in the file
   */
  async containerCount() {
    const e = await this.getSectionParsers(), { size: t } = await this.file.stat(), { cramContainerHeader1: r } = e;
    let a = 0, i = e.cramFileDefinition.maxLength;
    for (; i + r.maxLength + 8 < t; ) {
      const o = await this.getContainerAtPosition(i).getHeader();
      if (!o)
        break;
      if (a === 0) {
        i = o._endPosition;
        for (let s = 0; s < o.numBlocks; s += 1) {
          const c = await this.readBlock(i);
          if (c === void 0)
            return;
          i = c._endPosition;
        }
      } else
        i += o._size + o.length;
      a += 1;
    }
    return a;
  }
  getContainerAtPosition(e) {
    return new or(this, e);
  }
  async readBlockHeader(e) {
    const t = await this.getSectionParsers(), { cramBlockHeader: r } = t, { size: a } = await this.file.stat();
    if (e + r.maxLength >= a)
      return;
    const i = _.Buffer.allocUnsafe(r.maxLength);
    return await this.file.read(i, 0, r.maxLength, e), Y(i, r.parser, 0, e);
  }
  async _parseSection(e, t, r = e.maxLength, a = void 0) {
    let i;
    if (a)
      i = a;
    else {
      const { size: s } = await this.file.stat();
      if (t + r >= s)
        return;
      i = _.Buffer.allocUnsafe(r), await this.file.read(i, 0, r, t);
    }
    const o = Y(i, e.parser, 0, t);
    if (o._size !== r)
      throw new C(`section read error: requested size ${r} does not equal parsed size ${o._size}`);
    return o;
  }
  _uncompress(e, t, r) {
    if (e === "gzip")
      xt(t).copy(r);
    else if (e === "bzip2") {
      const a = bzip2.array(t);
      let i = bzip2.header(a), o = 0, s;
      do
        s = bzip2.decompress(a, i), s != -1 && (_.Buffer.from(s).copy(r, o), o += s.length, i -= s.length);
      while (s != -1);
    } else if (e === "rans")
      Nr(t, r);
    else if (e === "rans4x16")
      Ae.r4x16_uncompress(t, r);
    else if (e === "arith")
      Ae.arith_uncompress(t, r);
    else if (e === "fqzcomp")
      Ae.fqzcomp_uncompress(t, r);
    else if (e === "tok3")
      Ae.tok3_uncompress(t, r);
    else
      throw new G(`${e} decompression not yet implemented`);
  }
  async readBlock(e) {
    const { majorVersion: t } = await this.getDefinition(), r = await this.getSectionParsers(), a = await this.readBlockHeader(e);
    if (a === void 0)
      return;
    const i = a._endPosition, o = _.Buffer.allocUnsafe(a.uncompressedSize), s = {
      ...a,
      _endPosition: i,
      contentPosition: i,
      content: o
    };
    if (a.compressionMethod !== "raw") {
      const c = _.Buffer.allocUnsafe(a.compressedSize);
      await this.read(c, 0, a.compressedSize, i), this._uncompress(a.compressionMethod, c, o);
    } else
      await this.read(o, 0, a.uncompressedSize, i);
    if (t >= 3) {
      const c = await this._parseSection(r.cramBlockCrc32, i + a.compressedSize);
      if (c === void 0)
        return;
      s.crc32 = c.crc32, this.validateChecksums && await this.checkCrc32(e, a._size + a.compressedSize, c.crc32, "block data"), s._endPosition = c._endPosition, s._size = s.compressedSize + r.cramBlockCrc32.maxLength;
    } else
      s._endPosition = i + s.compressedSize, s._size = s.compressedSize;
    return s;
  }
}
"getDefinition getSectionParsers getSamHeader".split(" ").forEach((n) => Ie(rt, n));
class Ka {
  /**
   *
   * @param {object} args
   * @param {CramFile} args.cram
   * @param {Index-like} args.index object that supports getEntriesForRange(seqId,start,end) -> Promise[Array[index entries]]
   * @param {number} [args.cacheSize] optional maximum number of CRAM records to cache.  default 20,000
   * @param {number} [args.fetchSizeLimit] optional maximum number of bytes to fetch in a single getRecordsForRange call.  Default 3 MiB.
   * @param {boolean} [args.checkSequenceMD5] - default true. if false, disables verifying the MD5
   * checksum of the reference sequence underlying a slice. In some applications, this check can cause an inconvenient amount (many megabases) of sequences to be fetched.
   */
  constructor(e) {
    if (e.cram ? this.cram = e.cram : this.cram = new rt({
      url: e.cramUrl,
      path: e.cramPath,
      filehandle: e.cramFilehandle,
      seqFetch: e.seqFetch,
      checkSequenceMD5: e.checkSequenceMD5,
      cacheSize: e.cacheSize
    }), !(this.cram instanceof rt))
      throw new Error("invalid arguments: no cramfile");
    if (this.index = e.index, !this.index.getEntriesForRange)
      throw new Error("invalid arguments: not an index");
    this.fetchSizeLimit = e.fetchSizeLimit || 3e6;
  }
  /**
   *
   * @param {number} seq numeric ID of the reference sequence
   * @param {number} start start of the range of interest. 1-based closed coordinates.
   * @param {number} end end of the range of interest. 1-based closed coordinates.
   * @returns {Promise[Array[CramRecord]]}
   */
  async getRecordsForRange(e, t, r, a = {}) {
    if (a.viewAsPairs = a.viewAsPairs || !1, a.pairAcrossChr = a.pairAcrossChr || !1, a.maxInsertSize = a.maxInsertSize || 2e5, typeof e == "string")
      throw new G("string sequence names not yet supported");
    const i = e, o = await this.index.getEntriesForRange(i, t, r), s = o.map((d) => d.sliceBytes).reduce((d, h) => d + h, 0);
    if (s > this.fetchSizeLimit)
      throw new _r(`data size of ${s.toLocaleString()} bytes exceeded fetch size limit of ${this.fetchSizeLimit.toLocaleString()} bytes`);
    const c = (d) => d.sequenceId === e && d.alignmentStart <= r && d.lengthOnRef !== void 0 && d.alignmentStart + d.lengthOnRef - 1 >= t, l = await Promise.all(o.map((d) => this.getRecordsInSlice(d, c)));
    let f = Array.prototype.concat(...l);
    if (a.viewAsPairs) {
      const d = {}, h = {};
      for (let y = 0; y < f.length; y += 1) {
        const p = f[y].readName;
        if (p === void 0)
          throw new Error();
        const B = f[y].uniqueId;
        d[p] || (d[p] = 0), d[p] += 1, h[B] = 1;
      }
      const u = {};
      Object.entries(d).forEach(([y, p]) => {
        p === 1 && (u[y] = !0);
      });
      const g = [];
      for (let y = 0; y < f.length; y += 1) {
        const p = f[y], B = p.readName;
        if (B === void 0)
          throw new Error();
        if (u[B] && p.mate && (p.mate.sequenceId === i || a.pairAcrossChr) && Math.abs(p.alignmentStart - p.mate.alignmentStart) < a.maxInsertSize) {
          const k = this.index.getEntriesForRange(p.mate.sequenceId, p.mate.alignmentStart, p.mate.alignmentStart + 1);
          g.push(k);
        }
      }
      const m = await Promise.all(g);
      let v = [];
      for (let y = 0; y < m.length; y += 1)
        v.push(...m[y]);
      v = v.sort((y, p) => y.toString().localeCompare(p.toString())).filter((y, p, B) => !p || y.toString() !== B[p - 1].toString());
      const w = [], A = v.map((y) => y.sliceBytes).reduce((y, p) => y + p, 0);
      if (A > this.fetchSizeLimit)
        throw new Error(`mate data size of ${A.toLocaleString()} bytes exceeded fetch size limit of ${this.fetchSizeLimit.toLocaleString()} bytes`);
      v.forEach((y) => {
        let p = this.cram.featureCache.get(y.toString());
        p || (p = this.getRecordsInSlice(y, () => !0), this.cram.featureCache.set(y.toString(), p));
        const B = p.then((k) => {
          const I = [];
          for (let L = 0; L < k.length; L += 1) {
            const R = k[L];
            if (R.readName === void 0)
              throw new Error();
            u[R.readName] && !h[R.uniqueId] && I.push(R);
          }
          return I;
        });
        w.push(B);
      });
      const M = await Promise.all(w);
      if (M.length) {
        const y = M.reduce((p, B) => p.concat(B));
        f = f.concat(y);
      }
    }
    return f;
  }
  getRecordsInSlice({ containerStart: e, sliceStart: t, sliceBytes: r }, a) {
    return this.cram.getContainerAtPosition(e).getSlice(t, r).getRecords(a);
  }
  /**
   *
   * @param {number} seqId
   * @returns {Promise} true if the CRAM file contains data for the given
   * reference sequence numerical ID
   */
  hasDataForReferenceSequence(e) {
    return this.index.hasDataForReferenceSequence(e);
  }
}
const Xa = 21578050;
function Ct(n, e) {
  if (e.some((c) => c === void 0))
    throw new C("invalid .crai index file");
  const [t, r, a, i, o, s] = e;
  n[t] || (n[t] = []), n[t].push({
    start: r,
    span: a,
    containerStart: i,
    sliceStart: o,
    sliceBytes: s
  });
}
class Ya {
  /**
   *
   * @param {object} args
   * @param {string} [args.path]
   * @param {string} [args.url]
   * @param {FileHandle} [args.filehandle]
   */
  constructor(e) {
    this.filehandle = sr(e.url, e.path, e.filehandle), this._parseCache = new hr({
      cache: new Et({ maxSize: 1 }),
      fill: (t, r) => this.parseIndex()
    });
  }
  parseIndex() {
    const e = {};
    return this.filehandle.readFile().then((t) => t[0] === 31 && t[1] === 139 ? xt(t) : t).then((t) => {
      if (t.length > 4 && t.readUInt32LE(0) === Xa)
        throw new C("invalid .crai index file. note: file appears to be a .bai index. this is technically legal but please open a github issue if you need support");
      let r = [], a = "";
      for (let i = 0; i < t.length; i += 1) {
        const o = t[i];
        if (o >= 48 && o <= 57 || !a && o === 45)
          a += String.fromCharCode(o);
        else if (o === 9)
          r.push(Number.parseInt(a, 10)), a = "";
        else if (o === 10)
          r.push(Number.parseInt(a, 10)), a = "", Ct(e, r), r = [];
        else if (o !== 13 && o !== 32)
          throw new C("invalid .crai index file");
      }
      return a && r.push(Number.parseInt(a, 10)), r.length === 6 && Ct(e, r), Object.entries(e).forEach(([i, o]) => {
        e[i] = o.sort((s, c) => s.start - c.start || s.span - c.span);
      }), e;
    });
  }
  getIndex(e = {}) {
    return this._parseCache.get("index", null, e.signal);
  }
  /**
   * @param {number} seqId
   * @returns {Promise} true if the index contains entries for
   * the given reference sequence ID, false otherwise
   */
  async hasDataForReferenceSequence(e) {
    return !!(await this.getIndex())[e];
  }
  /**
   * fetch index entries for the given range
   *
   * @param {number} seqId
   * @param {number} queryStart
   * @param {number} queryEnd
   *
   * @returns {Promise} promise for
   * an array of objects of the form
   * `{start, span, containerStart, sliceStart, sliceBytes }`
   */
  async getEntriesForRange(e, t, r) {
    const a = (await this.getIndex())[e];
    if (!a)
      return [];
    const i = (s) => {
      const c = s.start, l = s.start + s.span;
      return c > r ? -1 : l <= t ? 1 : 0;
    }, o = [];
    for (let s = 0; s < a.length; s += 1)
      i(a[s]) === 0 && o.push(a[s]);
    return o;
  }
}
function Za(n, e, t) {
  if (!n)
    return [];
  const r = new Array(n.length);
  let a = 0, i = 0, o = 0, s = 0, c = e;
  for (let l = 0; l < n.length; l++) {
    const f = n[l], { code: d, pos: h, data: u, sub: g, ref: m } = f;
    if (s = o - c, c = o, s && i > 0 && (r[a++] = {
      start: o,
      type: "insertion",
      base: `${i}`,
      length: 0
    }, i = 0), o = f.refPos - 1 - e, d === "X")
      r[a++] = {
        start: o,
        length: 1,
        base: g,
        qual: t == null ? void 0 : t[h - 1],
        altbase: m == null ? void 0 : m.toUpperCase(),
        type: "mismatch"
      };
    else if (d === "I")
      r[a++] = {
        start: o,
        type: "insertion",
        base: `${u.length}`,
        length: 0
      };
    else if (d === "N")
      r[a++] = {
        type: "skip",
        length: u,
        start: o,
        base: "N"
      };
    else if (d === "S") {
      const v = u.length;
      r[a++] = {
        start: o,
        type: "softclip",
        base: `S${v}`,
        cliplen: v,
        length: 1
      };
    } else if (d !== "P")
      if (d === "H") {
        const v = u;
        r[a++] = {
          start: o,
          type: "hardclip",
          base: `H${v}`,
          cliplen: v,
          length: 1
        };
      } else
        d === "D" ? r[a++] = {
          type: "deletion",
          length: u,
          start: o,
          base: "*"
        } : d === "b" || d === "q" || d === "B" || d === "i" && i++;
  }
  return s && i > 0 && (r[a++] = {
    start: o,
    type: "insertion",
    base: `${i}`,
    length: 0
  }, i = 0), r.slice(0, a);
}
function Ja(n, e, t, r) {
  let a = "", i = "", o = "M", s = 0;
  if (!r)
    return "";
  const c = r.seq, l = r.start;
  let f = e, d = 0, h = 0;
  if (n !== void 0)
    for (let u = 0; u < n.length; u++) {
      const { code: g, refPos: m, sub: v, data: w } = n[u];
      if (d = m - f, a += c.slice(f - l, m - l), f = m, h > 0 && d && (i += `${h}I`, h = 0), s && o !== "M" && (i += `${s}${o}`, s = 0), d && (o = "M", s += d), g === "b") {
        const A = w.split(","), M = String.fromCharCode(...A);
        a += M, f += M.length, s += M.length;
      } else
        g === "B" || g === "X" ? (a += v, f++, s++) : g === "D" || g === "N" ? (f += w, s && (i += `${s}${o}`), i += w + g, s = 0) : g === "I" || g === "S" ? (a += w, s && (i += `${s}${o}`), i += w.length + g, s = 0) : g === "i" ? (s && (i += `${s}${o}`), h++, a += w, s = 0) : g === "P" ? (s && (i += `${s}${o}`), i += `${w}P`) : g === "H" && (s && (i += `${s}${o}`), i += `${w}H`, s = 0);
    }
  else
    d = t - a.length;
  return a.length !== t && (d = t - a.length, a += c.slice(f - l, f - l + d), s && o !== "M" && (i += `${s}${o}`, s = 0), o = "M", s += d), d && h > 0 && (i += `${h}I`), s && (i += `${s}${o}`), i;
}
class lt {
  // uses parameter properties to automatically create fields on the class
  // https://www.typescriptlang.org/docs/handbook/classes.html#parameter-properties
  constructor(e, t) {
    this.record = e, this._store = t;
  }
  _get_name() {
    return this.record.readName;
  }
  _get_start() {
    return this.record.alignmentStart - 1;
  }
  _get_end() {
    var e;
    return this.record.alignmentStart + ((e = this.record.lengthOnRef) !== null && e !== void 0 ? e : 1) - 1;
  }
  _get_cram_read_features() {
    return this.record.readFeatures;
  }
  _get_type() {
    return "match";
  }
  _get_score() {
    return this.record.mappingQuality;
  }
  _get_flags() {
    return this.record.flags;
  }
  _get_strand() {
    return this.record.isReverseComplemented() ? -1 : 1;
  }
  _read_group_id() {
    const e = this._store.samHeader.readGroups;
    return e ? e[this.record.readGroupId] : void 0;
  }
  _get_qual() {
    return (this.record.qualityScores || []).join(" ");
  }
  qualRaw() {
    return this.record.qualityScores;
  }
  _get_refName() {
    return this._store.refIdToName(this.record.sequenceId);
  }
  _get_is_paired() {
    return !!this.record.mate;
  }
  _get_pair_orientation() {
    return this.record.isPaired() ? this.record.getPairOrientation() : void 0;
  }
  _get_template_length() {
    return this.record.templateLength || this.record.templateSize;
  }
  _get_next_ref() {
    return this.record.mate ? this._store.refIdToName(this.record.mate.sequenceId) : void 0;
  }
  _get_next_segment_position() {
    return this.record.mate ? `${this._store.refIdToName(this.record.mate.sequenceId)}:${this.record.mate.alignmentStart}` : void 0;
  }
  _get_next_pos() {
    var e;
    return (e = this.record.mate) === null || e === void 0 ? void 0 : e.alignmentStart;
  }
  _get_tags() {
    const e = this._read_group_id(), { tags: t } = this.record;
    return e !== void 0 ? { ...t, RG: e } : t;
  }
  _get_seq() {
    return this.record.getReadBases();
  }
  // generate a CIGAR, based on code from jkbonfield
  _get_CIGAR() {
    return Ja(this.record.readFeatures, this.record.alignmentStart, this.record.readLength, this.record._refRegion);
  }
  tags() {
    return Object.getOwnPropertyNames(lt.prototype).filter((e) => e.startsWith("_get_") && e !== "_get_mismatches" && e !== "_get_cram_read_features").map((e) => e.replace("_get_", ""));
  }
  id() {
    return `${this._store.id}-${this.record.uniqueId}`;
  }
  get(e) {
    const t = `_get_${e}`;
    if (this[t])
      return this[t]();
  }
  parent() {
  }
  children() {
  }
  set() {
  }
  pairedFeature() {
    return !1;
  }
  _get_clipPos() {
    const e = this.get("mismatches");
    if (e.length) {
      const t = this.get("strand") === -1 ? e[e.length - 1] : e[0], { type: r, cliplen: a } = t;
      if (r === "softclip" || r === "hardclip")
        return a;
    }
    return 0;
  }
  toJSON() {
    return {
      ...Object.fromEntries(this.tags().map((e) => [e, this.get(e)]).filter((e) => e[1] !== void 0)),
      uniqueId: this.id()
    };
  }
  _get_mismatches() {
    const e = this.record.readFeatures, t = this.qualRaw(), r = this.get("start");
    return Za(e, r, t);
  }
}
class hi extends ur.BaseFeatureDataAdapter {
  constructor() {
    super(...arguments), this.samHeader = {}, this.seqIdToOriginalRefName = [];
  }
  async configurePre() {
    const e = this.getConf("cramLocation"), t = this.getConf("craiLocation");
    if (!e)
      throw new Error("missing cramLocation argument");
    if (!t)
      throw new Error("missing craiLocation argument");
    const r = this.pluginManager, a = new Ka({
      cramFilehandle: ht.openLocation(e, r),
      index: new Ya({ filehandle: ht.openLocation(t, r) }),
      seqFetch: (...s) => this.seqFetch(...s),
      checkSequenceMD5: !1,
      fetchSizeLimit: 2e8
      // just make this a large size to avoid hitting it
    });
    if (!this.getSubAdapter)
      throw new Error("Error getting subadapter");
    const i = this.getConf("sequenceAdapter"), o = await this.getSubAdapter(i);
    return {
      cram: a,
      sequenceAdapter: o.dataAdapter
    };
  }
  async configure() {
    return this.configureP || (this.configureP = this.configurePre().catch((e) => {
      throw this.configureP = void 0, e;
    })), this.configureP;
  }
  async getHeader(e) {
    const { cram: t } = await this.configure();
    return t.cram.getHeaderText();
  }
  async seqFetch(e, t, r) {
    t -= 1;
    const { sequenceAdapter: a } = await this.configure(), i = this.refIdToOriginalName(e) || this.refIdToName(e);
    if (!i)
      throw new Error("unknown");
    const s = (await gr.firstValueFrom(a.getFeatures({
      refName: i,
      start: t,
      end: r,
      assemblyName: ""
    }).pipe(vr.toArray()))).sort((c, l) => c.get("start") - l.get("start")).map((c) => {
      const l = c.get("start"), f = c.get("end"), d = Math.max(t - l, 0), u = Math.min(r - l, f - l) - d;
      return (c.get("seq") || c.get("residues")).slice(d, d + u);
    }).join("");
    if (s.length !== r - t)
      throw new Error(`sequence fetch failed: fetching ${i}:${(t - 1).toLocaleString()}-${r.toLocaleString()} returned ${s.length.toLocaleString()} bases, but should have returned ${(r - t).toLocaleString()}`);
    return s;
  }
  async setupPre(e) {
    const { statusCallback: t = () => {
    } } = e || {}, r = await this.configure();
    t("Downloading index");
    const { cram: a } = r, i = await a.cram.getSamHeader(), o = [], s = {};
    i.filter((f) => f.tag === "SQ").forEach((f, d) => {
      f.data.forEach((h) => {
        if (h.tag === "SN") {
          const u = h.value;
          s[u] = d, o[d] = u;
        }
      });
    });
    const c = i.filter((f) => f.tag === "RG").map((f) => {
      var d;
      return (d = f.data.find((h) => h.tag === "ID")) === null || d === void 0 ? void 0 : d.value;
    }), l = { idToName: o, nameToId: s, readGroups: c };
    return t(""), this.samHeader = l, { samHeader: l, ...r };
  }
  async setup(e) {
    return this.setupP || (this.setupP = this.setupPre(e).catch((t) => {
      throw this.setupP = void 0, t;
    })), this.setupP;
  }
  async getRefNames(e) {
    const { samHeader: t } = await this.setup(e);
    if (!t.idToName)
      throw new Error("CRAM file has no header lines");
    return t.idToName;
  }
  // use info from the SAM header if possible, but fall back to using
  // the ref seq order from when the browser's refseqs were loaded
  refNameToId(e) {
    if (this.samHeader.nameToId)
      return this.samHeader.nameToId[e];
    if (this.seqIdToRefName)
      return this.seqIdToRefName.indexOf(e);
  }
  // use info from the SAM header if possible, but fall back to using
  // the ref seq order from when the browser's refseqs were loaded
  refIdToName(e) {
    var t, r;
    return ((t = this.samHeader.idToName) === null || t === void 0 ? void 0 : t[e]) || ((r = this.seqIdToRefName) === null || r === void 0 ? void 0 : r[e]);
  }
  refIdToOriginalName(e) {
    return this.seqIdToOriginalRefName[e];
  }
  getFeatures(e, t) {
    const { signal: r, filterBy: a, statusCallback: i = () => {
    } } = t || {}, { refName: o, start: s, end: c, originalRefName: l } = e;
    return mr(async (f) => {
      const { cram: d } = await this.setup(t), h = this.refNameToId(o);
      if (h === void 0) {
        console.warn("Unknown refName", o), f.complete();
        return;
      }
      l && (this.seqIdToOriginalRefName[h] = l), i("Downloading alignments");
      const u = await d.getRecordsForRange(h, s, c);
      wr.checkAbortSignal(r);
      const { flagInclude: g = 0, flagExclude: m = 0, tagFilter: v, readName: w } = a || {};
      let A = u.filter((M) => {
        const y = M.flags;
        return (y & g) === g && !(y & m);
      });
      v && (A = A.filter((M) => {
        const y = M[v.tag];
        return y === "*" ? y !== void 0 : y === v.value;
      })), w && (A = A.filter((M) => M.readName === w)), A.forEach((M) => {
        f.next(this.cramRecordToFeature(M));
      }), i(""), f.complete();
    }, r);
  }
  freeResources() {
  }
  cramRecordToFeature(e) {
    return new lt(e, this);
  }
  // we return the configured fetchSizeLimit, and the bytes for the region
  async getMultiRegionFeatureDensityStats(e, t) {
    const r = await this.bytesForRegions(e, t), a = this.getConf("fetchSizeLimit");
    return {
      bytes: r,
      fetchSizeLimit: a
    };
  }
  /**
   * get the approximate number of bytes queried from the file for the given
   * query regions
   * @param regions - list of query regions
   */
  async bytesForRegions(e, t) {
    const { cram: r } = await this.configure();
    return (await Promise.all(e.map((i) => {
      const { refName: o, start: s, end: c } = i, l = this.refNameToId(o);
      return l !== void 0 ? r.index.getEntriesForRange(l, s, c) : [{ sliceBytes: 0 }];
    }))).flat().reduce((i, o) => i + o.sliceBytes, 0);
  }
}
export {
  hi as default
};
