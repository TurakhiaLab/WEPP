import { b8 as rn, b3 as Ke, b6 as sn, b5 as on, aX as Rt, aY as fn } from "./index-8ee3baa4.js";
import { Q as ln } from "./index-05a6222b.js";
const qe = BigInt(32);
function hn(t, i, e) {
  const n = +!!e, a = +!e;
  return BigInt(t.getInt32(i, e) * a + t.getInt32(i + 4, e) * n) << qe | BigInt(t.getUint32(i, e) * n + t.getUint32(i + 4, e) * a);
}
function un(t, i, e) {
  const n = t.getUint32(i, e), a = t.getUint32(i + 4, e), r = +!!e, s = +!e;
  return BigInt(n * s + a * r) << qe | BigInt(n * r + a * s);
}
function dn(t, i, e, n) {
  const a = Number(e >> qe), r = Number(e & BigInt(4294967295));
  n ? (t.setInt32(i + 4, a, n), t.setUint32(i, r, n)) : (t.setInt32(i, a, n), t.setUint32(i + 4, r, n));
}
function cn(t, i, e, n) {
  const a = Number(e >> qe), r = Number(e & BigInt(4294967295));
  n ? (t.setUint32(i + 4, a, n), t.setUint32(i, r, n)) : (t.setUint32(i, a, n), t.setUint32(i + 4, r, n));
}
"getBigInt64" in DataView || (DataView.prototype.getBigInt64 = function(t, i) {
  return hn(this, t, i);
});
"getBigUint64" in DataView || (DataView.prototype.getBigUint64 = function(t, i) {
  return un(this, t, i);
});
"setBigInt64" in DataView || (DataView.prototype.setBigInt64 = function(t, i, e) {
  dn(this, t, i, e);
});
"setBigUint64" in DataView || (DataView.prototype.setBigUint64 = function(t, i, e) {
  cn(this, t, i, e);
});
class _n {
  constructor(i, e) {
    this.code = "", this.scopes = [["vars"]], this.bitFields = [], this.tmpVariableCount = 0, this.references = /* @__PURE__ */ new Map(), this.imports = [], this.reverseImports = /* @__PURE__ */ new Map(), this.useContextVariables = !1, this.importPath = i, this.useContextVariables = e;
  }
  generateVariable(i) {
    const e = [...this.scopes[this.scopes.length - 1]];
    return i && e.push(i), e.join(".");
  }
  generateOption(i) {
    switch (typeof i) {
      case "number":
        return i.toString();
      case "string":
        return this.generateVariable(i);
      case "function":
        return `${this.addImport(i)}.call(${this.generateVariable()}, vars)`;
    }
  }
  generateError(i) {
    this.pushCode(`throw new Error(${i});`);
  }
  generateTmpVariable() {
    return "$tmp" + this.tmpVariableCount++;
  }
  pushCode(i) {
    this.code += i + `
`;
  }
  pushPath(i) {
    i && this.scopes[this.scopes.length - 1].push(i);
  }
  popPath(i) {
    i && this.scopes[this.scopes.length - 1].pop();
  }
  pushScope(i) {
    this.scopes.push([i]);
  }
  popScope() {
    this.scopes.pop();
  }
  addImport(i) {
    if (!this.importPath)
      return `(${i})`;
    let e = this.reverseImports.get(i);
    return e || (e = this.imports.push(i) - 1, this.reverseImports.set(i, e)), `${this.importPath}[${e}]`;
  }
  addReference(i) {
    this.references.has(i) || this.references.set(i, { resolved: !1, requested: !1 });
  }
  markResolved(i) {
    const e = this.references.get(i);
    e && (e.resolved = !0);
  }
  markRequested(i) {
    i.forEach((e) => {
      const n = this.references.get(e);
      n && (n.requested = !0);
    });
  }
  getUnresolvedReferences() {
    return Array.from(this.references).filter(([i, e]) => !e.resolved && !e.requested).map(([i, e]) => i);
  }
}
const Z = /* @__PURE__ */ new Map(), ie = "___parser_", L = {
  uint8: 1,
  uint16le: 2,
  uint16be: 2,
  uint32le: 4,
  uint32be: 4,
  int8: 1,
  int16le: 2,
  int16be: 2,
  int32le: 4,
  int32be: 4,
  int64be: 8,
  int64le: 8,
  uint64be: 8,
  uint64le: 8,
  floatle: 4,
  floatbe: 4,
  doublele: 8,
  doublebe: 8
}, Ve = {
  uint8: "Uint8",
  uint16le: "Uint16",
  uint16be: "Uint16",
  uint32le: "Uint32",
  uint32be: "Uint32",
  int8: "Int8",
  int16le: "Int16",
  int16be: "Int16",
  int32le: "Int32",
  int32be: "Int32",
  int64be: "BigInt64",
  int64le: "BigInt64",
  uint64be: "BigUint64",
  uint64le: "BigUint64",
  floatle: "Float32",
  floatbe: "Float32",
  doublele: "Float64",
  doublebe: "Float64"
}, Fe = {
  uint8: !1,
  uint16le: !0,
  uint16be: !1,
  uint32le: !0,
  uint32be: !1,
  int8: !1,
  int16le: !0,
  int16be: !1,
  int32le: !0,
  int32be: !1,
  int64be: !1,
  int64le: !0,
  uint64be: !1,
  uint64le: !0,
  floatle: !0,
  floatbe: !1,
  doublele: !0,
  doublebe: !1
};
class x {
  constructor() {
    this.varName = "", this.type = "", this.options = {}, this.endian = "be", this.useContextVariables = !1;
  }
  static start() {
    return new x();
  }
  primitiveGenerateN(i, e) {
    const n = Ve[i], a = Fe[i];
    e.pushCode(`${e.generateVariable(this.varName)} = dataView.get${n}(offset, ${a});`), e.pushCode(`offset += ${L[i]};`);
  }
  primitiveN(i, e, n) {
    return this.setNextParser(i, e, n);
  }
  useThisEndian(i) {
    return i + this.endian.toLowerCase();
  }
  uint8(i, e = {}) {
    return this.primitiveN("uint8", i, e);
  }
  uint16(i, e = {}) {
    return this.primitiveN(this.useThisEndian("uint16"), i, e);
  }
  uint16le(i, e = {}) {
    return this.primitiveN("uint16le", i, e);
  }
  uint16be(i, e = {}) {
    return this.primitiveN("uint16be", i, e);
  }
  uint32(i, e = {}) {
    return this.primitiveN(this.useThisEndian("uint32"), i, e);
  }
  uint32le(i, e = {}) {
    return this.primitiveN("uint32le", i, e);
  }
  uint32be(i, e = {}) {
    return this.primitiveN("uint32be", i, e);
  }
  int8(i, e = {}) {
    return this.primitiveN("int8", i, e);
  }
  int16(i, e = {}) {
    return this.primitiveN(this.useThisEndian("int16"), i, e);
  }
  int16le(i, e = {}) {
    return this.primitiveN("int16le", i, e);
  }
  int16be(i, e = {}) {
    return this.primitiveN("int16be", i, e);
  }
  int32(i, e = {}) {
    return this.primitiveN(this.useThisEndian("int32"), i, e);
  }
  int32le(i, e = {}) {
    return this.primitiveN("int32le", i, e);
  }
  int32be(i, e = {}) {
    return this.primitiveN("int32be", i, e);
  }
  bigIntVersionCheck() {
    if (!DataView.prototype.getBigInt64)
      throw new Error("BigInt64 is unsupported on this runtime");
  }
  int64(i, e = {}) {
    return this.bigIntVersionCheck(), this.primitiveN(this.useThisEndian("int64"), i, e);
  }
  int64be(i, e = {}) {
    return this.bigIntVersionCheck(), this.primitiveN("int64be", i, e);
  }
  int64le(i, e = {}) {
    return this.bigIntVersionCheck(), this.primitiveN("int64le", i, e);
  }
  uint64(i, e = {}) {
    return this.bigIntVersionCheck(), this.primitiveN(this.useThisEndian("uint64"), i, e);
  }
  uint64be(i, e = {}) {
    return this.bigIntVersionCheck(), this.primitiveN("uint64be", i, e);
  }
  uint64le(i, e = {}) {
    return this.bigIntVersionCheck(), this.primitiveN("uint64le", i, e);
  }
  floatle(i, e = {}) {
    return this.primitiveN("floatle", i, e);
  }
  floatbe(i, e = {}) {
    return this.primitiveN("floatbe", i, e);
  }
  doublele(i, e = {}) {
    return this.primitiveN("doublele", i, e);
  }
  doublebe(i, e = {}) {
    return this.primitiveN("doublebe", i, e);
  }
  bitN(i, e, n) {
    return n.length = i, this.setNextParser("bit", e, n);
  }
  bit1(i, e = {}) {
    return this.bitN(1, i, e);
  }
  bit2(i, e = {}) {
    return this.bitN(2, i, e);
  }
  bit3(i, e = {}) {
    return this.bitN(3, i, e);
  }
  bit4(i, e = {}) {
    return this.bitN(4, i, e);
  }
  bit5(i, e = {}) {
    return this.bitN(5, i, e);
  }
  bit6(i, e = {}) {
    return this.bitN(6, i, e);
  }
  bit7(i, e = {}) {
    return this.bitN(7, i, e);
  }
  bit8(i, e = {}) {
    return this.bitN(8, i, e);
  }
  bit9(i, e = {}) {
    return this.bitN(9, i, e);
  }
  bit10(i, e = {}) {
    return this.bitN(10, i, e);
  }
  bit11(i, e = {}) {
    return this.bitN(11, i, e);
  }
  bit12(i, e = {}) {
    return this.bitN(12, i, e);
  }
  bit13(i, e = {}) {
    return this.bitN(13, i, e);
  }
  bit14(i, e = {}) {
    return this.bitN(14, i, e);
  }
  bit15(i, e = {}) {
    return this.bitN(15, i, e);
  }
  bit16(i, e = {}) {
    return this.bitN(16, i, e);
  }
  bit17(i, e = {}) {
    return this.bitN(17, i, e);
  }
  bit18(i, e = {}) {
    return this.bitN(18, i, e);
  }
  bit19(i, e = {}) {
    return this.bitN(19, i, e);
  }
  bit20(i, e = {}) {
    return this.bitN(20, i, e);
  }
  bit21(i, e = {}) {
    return this.bitN(21, i, e);
  }
  bit22(i, e = {}) {
    return this.bitN(22, i, e);
  }
  bit23(i, e = {}) {
    return this.bitN(23, i, e);
  }
  bit24(i, e = {}) {
    return this.bitN(24, i, e);
  }
  bit25(i, e = {}) {
    return this.bitN(25, i, e);
  }
  bit26(i, e = {}) {
    return this.bitN(26, i, e);
  }
  bit27(i, e = {}) {
    return this.bitN(27, i, e);
  }
  bit28(i, e = {}) {
    return this.bitN(28, i, e);
  }
  bit29(i, e = {}) {
    return this.bitN(29, i, e);
  }
  bit30(i, e = {}) {
    return this.bitN(30, i, e);
  }
  bit31(i, e = {}) {
    return this.bitN(31, i, e);
  }
  bit32(i, e = {}) {
    return this.bitN(32, i, e);
  }
  namely(i) {
    return Z.set(i, this), this.alias = i, this;
  }
  skip(i, e = {}) {
    return this.seek(i, e);
  }
  seek(i, e = {}) {
    if (e.assert)
      throw new Error("assert option on seek is not allowed.");
    return this.setNextParser("seek", "", { length: i });
  }
  string(i, e) {
    if (!e.zeroTerminated && !e.length && !e.greedy)
      throw new Error("One of length, zeroTerminated, or greedy must be defined for string.");
    if ((e.zeroTerminated || e.length) && e.greedy)
      throw new Error("greedy is mutually exclusive with length and zeroTerminated for string.");
    if (e.stripNull && !(e.length || e.greedy))
      throw new Error("length or greedy must be defined if stripNull is enabled.");
    return e.encoding = e.encoding || "utf8", this.setNextParser("string", i, e);
  }
  buffer(i, e) {
    if (!e.length && !e.readUntil)
      throw new Error("length or readUntil must be defined for buffer.");
    return this.setNextParser("buffer", i, e);
  }
  wrapped(i, e) {
    if (typeof e != "object" && typeof i == "object" && (e = i, i = ""), !e || !e.wrapper || !e.type)
      throw new Error("Both wrapper and type must be defined for wrapped.");
    if (!e.length && !e.readUntil)
      throw new Error("length or readUntil must be defined for wrapped.");
    return this.setNextParser("wrapper", i, e);
  }
  array(i, e) {
    if (!e.readUntil && !e.length && !e.lengthInBytes)
      throw new Error("One of readUntil, length and lengthInBytes must be defined for array.");
    if (!e.type)
      throw new Error("type is required for array.");
    if (typeof e.type == "string" && !Z.has(e.type) && !(e.type in L))
      throw new Error(`Array element type "${e.type}" is unkown.`);
    return this.setNextParser("array", i, e);
  }
  choice(i, e) {
    if (typeof e != "object" && typeof i == "object" && (e = i, i = ""), !e)
      throw new Error("tag and choices are are required for choice.");
    if (!e.tag)
      throw new Error("tag is requird for choice.");
    if (!e.choices)
      throw new Error("choices is required for choice.");
    for (const n in e.choices) {
      const a = parseInt(n, 10), r = e.choices[a];
      if (isNaN(a))
        throw new Error(`Choice key "${n}" is not a number.`);
      if (typeof r == "string" && !Z.has(r) && !(r in L))
        throw new Error(`Choice type "${r}" is unkown.`);
    }
    return this.setNextParser("choice", i, e);
  }
  nest(i, e) {
    if (typeof e != "object" && typeof i == "object" && (e = i, i = ""), !e || !e.type)
      throw new Error("type is required for nest.");
    if (!(e.type instanceof x) && !Z.has(e.type))
      throw new Error("type must be a known parser name or a Parser object.");
    if (!(e.type instanceof x) && !i)
      throw new Error("type must be a Parser object if the variable name is omitted.");
    return this.setNextParser("nest", i, e);
  }
  pointer(i, e) {
    if (!e.offset)
      throw new Error("offset is required for pointer.");
    if (!e.type)
      throw new Error("type is required for pointer.");
    if (typeof e.type == "string" && !(e.type in L) && !Z.has(e.type))
      throw new Error(`Pointer type "${e.type}" is unkown.`);
    return this.setNextParser("pointer", i, e);
  }
  saveOffset(i, e = {}) {
    return this.setNextParser("saveOffset", i, e);
  }
  endianness(i) {
    switch (i.toLowerCase()) {
      case "little":
        this.endian = "le";
        break;
      case "big":
        this.endian = "be";
        break;
      default:
        throw new Error('endianness must be one of "little" or "big"');
    }
    return this;
  }
  endianess(i) {
    return this.endianness(i);
  }
  useContextVars(i = !0) {
    return this.useContextVariables = i, this;
  }
  create(i) {
    if (!(i instanceof Function))
      throw new Error("Constructor must be a Function object.");
    return this.constructorFn = i, this;
  }
  getContext(i) {
    const e = new _n(i, this.useContextVariables);
    return e.pushCode("var dataView = new DataView(buffer.buffer, buffer.byteOffset, buffer.length);"), this.alias ? (this.addAliasedCode(e), e.pushCode(`return ${ie + this.alias}(0).result;`)) : this.addRawCode(e), e;
  }
  getCode() {
    const i = "imports";
    return this.getContext(i).code;
  }
  addRawCode(i) {
    i.pushCode("var offset = 0;"), i.pushCode(`var vars = ${this.constructorFn ? "new constructorFn()" : "{}"};`), i.pushCode("vars.$parent = null;"), i.pushCode("vars.$root = vars;"), this.generate(i), this.resolveReferences(i), i.pushCode("delete vars.$parent;"), i.pushCode("delete vars.$root;"), i.pushCode("return vars;");
  }
  addAliasedCode(i) {
    return i.pushCode(`function ${ie + this.alias}(offset, context) {`), i.pushCode(`var vars = ${this.constructorFn ? "new constructorFn()" : "{}"};`), i.pushCode("var ctx = Object.assign({$parent: null, $root: vars}, context || {});"), i.pushCode("vars = Object.assign(vars, ctx);"), this.generate(i), i.markResolved(this.alias), this.resolveReferences(i), i.pushCode("Object.keys(ctx).forEach(function (item) { delete vars[item]; });"), i.pushCode("return { offset: offset, result: vars };"), i.pushCode("}"), i;
  }
  resolveReferences(i) {
    const e = i.getUnresolvedReferences();
    i.markRequested(e), e.forEach((n) => {
      var a;
      (a = Z.get(n)) === null || a === void 0 || a.addAliasedCode(i);
    });
  }
  compile() {
    const i = "imports", e = this.getContext(i);
    this.compiled = new Function(i, "TextDecoder", `return function (buffer, constructorFn) { ${e.code} };`)(e.imports, TextDecoder);
  }
  sizeOf() {
    let i = NaN;
    if (Object.keys(L).indexOf(this.type) >= 0)
      i = L[this.type];
    else if (this.type === "string" && typeof this.options.length == "number")
      i = this.options.length;
    else if (this.type === "buffer" && typeof this.options.length == "number")
      i = this.options.length;
    else if (this.type === "array" && typeof this.options.length == "number") {
      let e = NaN;
      typeof this.options.type == "string" ? e = L[this.options.type] : this.options.type instanceof x && (e = this.options.type.sizeOf()), i = this.options.length * e;
    } else
      this.type === "seek" ? i = this.options.length : this.type === "nest" ? i = this.options.type.sizeOf() : this.type || (i = 0);
    return this.next && (i += this.next.sizeOf()), i;
  }
  // Follow the parser chain till the root and start parsing from there
  parse(i) {
    return this.compiled || this.compile(), this.compiled(i, this.constructorFn);
  }
  setNextParser(i, e, n) {
    const a = new x();
    return a.type = i, a.varName = e, a.options = n, a.endian = this.endian, this.head ? this.head.next = a : this.next = a, this.head = a, this;
  }
  // Call code generator for this parser
  generate(i) {
    if (this.type) {
      switch (this.type) {
        case "uint8":
        case "uint16le":
        case "uint16be":
        case "uint32le":
        case "uint32be":
        case "int8":
        case "int16le":
        case "int16be":
        case "int32le":
        case "int32be":
        case "int64be":
        case "int64le":
        case "uint64be":
        case "uint64le":
        case "floatle":
        case "floatbe":
        case "doublele":
        case "doublebe":
          this.primitiveGenerateN(this.type, i);
          break;
        case "bit":
          this.generateBit(i);
          break;
        case "string":
          this.generateString(i);
          break;
        case "buffer":
          this.generateBuffer(i);
          break;
        case "seek":
          this.generateSeek(i);
          break;
        case "nest":
          this.generateNest(i);
          break;
        case "array":
          this.generateArray(i);
          break;
        case "choice":
          this.generateChoice(i);
          break;
        case "pointer":
          this.generatePointer(i);
          break;
        case "saveOffset":
          this.generateSaveOffset(i);
          break;
        case "wrapper":
          this.generateWrapper(i);
          break;
      }
      this.type !== "bit" && this.generateAssert(i);
    }
    const e = i.generateVariable(this.varName);
    return this.options.formatter && this.type !== "bit" && this.generateFormatter(i, e, this.options.formatter), this.generateNext(i);
  }
  generateAssert(i) {
    if (!this.options.assert)
      return;
    const e = i.generateVariable(this.varName);
    switch (typeof this.options.assert) {
      case "function":
        {
          const n = i.addImport(this.options.assert);
          i.pushCode(`if (!${n}.call(vars, ${e})) {`);
        }
        break;
      case "number":
        i.pushCode(`if (${this.options.assert} !== ${e}) {`);
        break;
      case "string":
        i.pushCode(`if (${JSON.stringify(this.options.assert)} !== ${e}) {`);
        break;
      default:
        throw new Error("assert option must be a string, number or a function.");
    }
    i.generateError(`"Assertion error: ${e} is " + ${JSON.stringify(this.options.assert.toString())}`), i.pushCode("}");
  }
  // Recursively call code generators and append results
  generateNext(i) {
    return this.next && (i = this.next.generate(i)), i;
  }
  generateBit(i) {
    const e = JSON.parse(JSON.stringify(this));
    if (e.options = this.options, e.generateAssert = this.generateAssert.bind(this), e.generateFormatter = this.generateFormatter.bind(this), e.varName = i.generateVariable(e.varName), i.bitFields.push(e), !this.next || this.next && ["bit", "nest"].indexOf(this.next.type) < 0) {
      const n = i.generateTmpVariable();
      i.pushCode(`var ${n} = 0;`);
      const a = (l = 0) => {
        let b = 0;
        for (let d = l; d < i.bitFields.length; d++) {
          const h = i.bitFields[d].options.length;
          if (b + h > 32)
            break;
          b += h;
        }
        return b;
      }, r = (l) => (l <= 8 ? (i.pushCode(`${n} = dataView.getUint8(offset);`), l = 8) : l <= 16 ? (i.pushCode(`${n} = dataView.getUint16(offset);`), l = 16) : l <= 24 ? (i.pushCode(`${n} = (dataView.getUint16(offset) << 8) | dataView.getUint8(offset + 2);`), l = 24) : (i.pushCode(`${n} = dataView.getUint32(offset);`), l = 32), i.pushCode(`offset += ${l / 8};`), l);
      let s = 0;
      const o = this.endian === "be";
      let u = 0, f = 0;
      i.bitFields.forEach((l, b) => {
        let d = l.options.length;
        if (d > f) {
          if (f) {
            const k = -1 >>> 32 - f;
            i.pushCode(`${l.varName} = (${n} & 0x${k.toString(16)}) << ${d - f};`), d -= f;
          }
          s = 0, f = u = r(a(b) - f);
        }
        const h = o ? u - s - d : s, p = -1 >>> 32 - d;
        i.pushCode(`${l.varName} ${d < l.options.length ? "|=" : "="} ${n} >> ${h} & 0x${p.toString(16)};`), l.options.length === 32 && i.pushCode(`${l.varName} >>>= 0`), l.options.assert && l.generateAssert(i), l.options.formatter && l.generateFormatter(i, l.varName, l.options.formatter), s += d, f -= d;
      }), i.bitFields = [];
    }
  }
  generateSeek(i) {
    const e = i.generateOption(this.options.length);
    i.pushCode(`offset += ${e};`);
  }
  generateString(i) {
    const e = i.generateVariable(this.varName), n = i.generateTmpVariable(), a = this.options.encoding, r = a.toLowerCase() === "hex", s = 'b => b.toString(16).padStart(2, "0")';
    if (this.options.length && this.options.zeroTerminated) {
      const o = this.options.length;
      i.pushCode(`var ${n} = offset;`), i.pushCode(`while(dataView.getUint8(offset++) !== 0 && offset - ${n} < ${o});`);
      const u = `offset - ${n} < ${o} ? offset - 1 : offset`;
      i.pushCode(r ? `${e} = Array.from(buffer.subarray(${n}, ${u}), ${s}).join('');` : `${e} = new TextDecoder('${a}').decode(buffer.subarray(${n}, ${u}));`);
    } else if (this.options.length) {
      const o = i.generateOption(this.options.length);
      i.pushCode(r ? `${e} = Array.from(buffer.subarray(offset, offset + ${o}), ${s}).join('');` : `${e} = new TextDecoder('${a}').decode(buffer.subarray(offset, offset + ${o}));`), i.pushCode(`offset += ${o};`);
    } else
      this.options.zeroTerminated ? (i.pushCode(`var ${n} = offset;`), i.pushCode("while(dataView.getUint8(offset++) !== 0);"), i.pushCode(r ? `${e} = Array.from(buffer.subarray(${n}, offset - 1), ${s}).join('');` : `${e} = new TextDecoder('${a}').decode(buffer.subarray(${n}, offset - 1));`)) : this.options.greedy && (i.pushCode(`var ${n} = offset;`), i.pushCode("while(buffer.length > offset++);"), i.pushCode(r ? `${e} = Array.from(buffer.subarray(${n}, offset), ${s}).join('');` : `${e} = new TextDecoder('${a}').decode(buffer.subarray(${n}, offset));`));
    this.options.stripNull && i.pushCode(`${e} = ${e}.replace(/\\x00+$/g, '')`);
  }
  generateBuffer(i) {
    const e = i.generateVariable(this.varName);
    if (typeof this.options.readUntil == "function") {
      const n = this.options.readUntil, a = i.generateTmpVariable(), r = i.generateTmpVariable();
      i.pushCode(`var ${a} = offset;`), i.pushCode(`var ${r} = 0;`), i.pushCode("while (offset < buffer.length) {"), i.pushCode(`${r} = dataView.getUint8(offset);`);
      const s = i.addImport(n);
      i.pushCode(`if (${s}.call(${i.generateVariable()}, ${r}, buffer.subarray(offset))) break;`), i.pushCode("offset += 1;"), i.pushCode("}"), i.pushCode(`${e} = buffer.subarray(${a}, offset);`);
    } else if (this.options.readUntil === "eof")
      i.pushCode(`${e} = buffer.subarray(offset);`);
    else {
      const n = i.generateOption(this.options.length);
      i.pushCode(`${e} = buffer.subarray(offset, offset + ${n});`), i.pushCode(`offset += ${n};`);
    }
    this.options.clone && i.pushCode(`${e} = buffer.constructor.from(${e});`);
  }
  generateArray(i) {
    const e = i.generateOption(this.options.length), n = i.generateOption(this.options.lengthInBytes), a = this.options.type, r = i.generateTmpVariable(), s = i.generateVariable(this.varName), o = i.generateTmpVariable(), u = this.options.key, f = typeof u == "string";
    if (f ? i.pushCode(`${s} = {};`) : i.pushCode(`${s} = [];`), typeof this.options.readUntil == "function" ? i.pushCode("do {") : this.options.readUntil === "eof" ? i.pushCode(`for (var ${r} = 0; offset < buffer.length; ${r}++) {`) : n !== void 0 ? i.pushCode(`for (var ${r} = offset + ${n}; offset < ${r}; ) {`) : i.pushCode(`for (var ${r} = ${e}; ${r} > 0; ${r}--) {`), typeof a == "string")
      if (Z.get(a)) {
        const l = i.generateTmpVariable();
        if (i.pushCode(`var ${l} = ${ie + a}(offset, {`), i.useContextVariables) {
          const b = i.generateVariable();
          i.pushCode(`$parent: ${b},`), i.pushCode(`$root: ${b}.$root,`), !this.options.readUntil && n === void 0 && i.pushCode(`$index: ${e} - ${r},`);
        }
        i.pushCode("});"), i.pushCode(`var ${o} = ${l}.result; offset = ${l}.offset;`), a !== this.alias && i.addReference(a);
      } else {
        const l = Ve[a], b = Fe[a];
        i.pushCode(`var ${o} = dataView.get${l}(offset, ${b});`), i.pushCode(`offset += ${L[a]};`);
      }
    else if (a instanceof x) {
      i.pushCode(`var ${o} = {};`);
      const l = i.generateVariable();
      i.pushScope(o), i.useContextVariables && (i.pushCode(`${o}.$parent = ${l};`), i.pushCode(`${o}.$root = ${l}.$root;`), !this.options.readUntil && n === void 0 && i.pushCode(`${o}.$index = ${e} - ${r};`)), a.generate(i), i.useContextVariables && (i.pushCode(`delete ${o}.$parent;`), i.pushCode(`delete ${o}.$root;`), i.pushCode(`delete ${o}.$index;`)), i.popScope();
    }
    if (f ? i.pushCode(`${s}[${o}.${u}] = ${o};`) : i.pushCode(`${s}.push(${o});`), i.pushCode("}"), typeof this.options.readUntil == "function") {
      const l = this.options.readUntil, b = i.addImport(l);
      i.pushCode(`while (!${b}.call(${i.generateVariable()}, ${o}, buffer.subarray(offset)));`);
    }
  }
  generateChoiceCase(i, e, n) {
    if (typeof n == "string") {
      const a = i.generateVariable(this.varName);
      if (Z.has(n)) {
        const r = i.generateTmpVariable();
        i.pushCode(`var ${r} = ${ie + n}(offset, {`), i.useContextVariables && (i.pushCode(`$parent: ${a}.$parent,`), i.pushCode(`$root: ${a}.$root,`)), i.pushCode("});"), i.pushCode(`${a} = ${r}.result; offset = ${r}.offset;`), n !== this.alias && i.addReference(n);
      } else {
        const r = Ve[n], s = Fe[n];
        i.pushCode(`${a} = dataView.get${r}(offset, ${s});`), i.pushCode(`offset += ${L[n]}`);
      }
    } else
      n instanceof x && (i.pushPath(e), n.generate(i), i.popPath(e));
  }
  generateChoice(i) {
    const e = i.generateOption(this.options.tag), n = i.generateVariable(this.varName);
    if (this.varName && (i.pushCode(`${n} = {};`), i.useContextVariables)) {
      const a = i.generateVariable();
      i.pushCode(`${n}.$parent = ${a};`), i.pushCode(`${n}.$root = ${a}.$root;`);
    }
    i.pushCode(`switch(${e}) {`);
    for (const a in this.options.choices) {
      const r = parseInt(a, 10), s = this.options.choices[r];
      i.pushCode(`case ${r}:`), this.generateChoiceCase(i, this.varName, s), i.pushCode("break;");
    }
    i.pushCode("default:"), this.options.defaultChoice ? this.generateChoiceCase(i, this.varName, this.options.defaultChoice) : i.generateError(`"Met undefined tag value " + ${e} + " at choice"`), i.pushCode("}"), this.varName && i.useContextVariables && (i.pushCode(`delete ${n}.$parent;`), i.pushCode(`delete ${n}.$root;`));
  }
  generateNest(i) {
    const e = i.generateVariable(this.varName);
    if (this.options.type instanceof x) {
      if (this.varName && (i.pushCode(`${e} = {};`), i.useContextVariables)) {
        const n = i.generateVariable();
        i.pushCode(`${e}.$parent = ${n};`), i.pushCode(`${e}.$root = ${n}.$root;`);
      }
      i.pushPath(this.varName), this.options.type.generate(i), i.popPath(this.varName), this.varName && i.useContextVariables && i.useContextVariables && (i.pushCode(`delete ${e}.$parent;`), i.pushCode(`delete ${e}.$root;`));
    } else if (Z.has(this.options.type)) {
      const n = i.generateTmpVariable();
      if (i.pushCode(`var ${n} = ${ie + this.options.type}(offset, {`), i.useContextVariables) {
        const a = i.generateVariable();
        i.pushCode(`$parent: ${a},`), i.pushCode(`$root: ${a}.$root,`);
      }
      i.pushCode("});"), i.pushCode(`${e} = ${n}.result; offset = ${n}.offset;`), this.options.type !== this.alias && i.addReference(this.options.type);
    }
  }
  generateWrapper(i) {
    const e = i.generateVariable(this.varName), n = i.generateTmpVariable();
    if (typeof this.options.readUntil == "function") {
      const u = this.options.readUntil, f = i.generateTmpVariable(), l = i.generateTmpVariable();
      i.pushCode(`var ${f} = offset;`), i.pushCode(`var ${l} = 0;`), i.pushCode("while (offset < buffer.length) {"), i.pushCode(`${l} = dataView.getUint8(offset);`);
      const b = i.addImport(u);
      i.pushCode(`if (${b}.call(${i.generateVariable()}, ${l}, buffer.subarray(offset))) break;`), i.pushCode("offset += 1;"), i.pushCode("}"), i.pushCode(`${n} = buffer.subarray(${f}, offset);`);
    } else if (this.options.readUntil === "eof")
      i.pushCode(`${n} = buffer.subarray(offset);`);
    else {
      const u = i.generateOption(this.options.length);
      i.pushCode(`${n} = buffer.subarray(offset, offset + ${u});`), i.pushCode(`offset += ${u};`);
    }
    this.options.clone && i.pushCode(`${n} = buffer.constructor.from(${n});`);
    const a = i.generateTmpVariable(), r = i.generateTmpVariable(), s = i.generateTmpVariable(), o = i.addImport(this.options.wrapper);
    if (i.pushCode(`${n} = ${o}.call(this, ${n}).subarray(0);`), i.pushCode(`var ${a} = buffer;`), i.pushCode(`var ${r} = offset;`), i.pushCode(`var ${s} = dataView;`), i.pushCode(`buffer = ${n};`), i.pushCode("offset = 0;"), i.pushCode("dataView = new DataView(buffer.buffer, buffer.byteOffset, buffer.length);"), this.options.type instanceof x)
      this.varName && i.pushCode(`${e} = {};`), i.pushPath(this.varName), this.options.type.generate(i), i.popPath(this.varName);
    else if (Z.has(this.options.type)) {
      const u = i.generateTmpVariable();
      i.pushCode(`var ${u} = ${ie + this.options.type}(0);`), i.pushCode(`${e} = ${u}.result;`), this.options.type !== this.alias && i.addReference(this.options.type);
    }
    i.pushCode(`buffer = ${a};`), i.pushCode(`dataView = ${s};`), i.pushCode(`offset = ${r};`);
  }
  generateFormatter(i, e, n) {
    if (typeof n == "function") {
      const a = i.addImport(n);
      i.pushCode(`${e} = ${a}.call(${i.generateVariable()}, ${e});`);
    }
  }
  generatePointer(i) {
    const e = this.options.type, n = i.generateOption(this.options.offset), a = i.generateTmpVariable(), r = i.generateVariable(this.varName);
    if (i.pushCode(`var ${a} = offset;`), i.pushCode(`offset = ${n};`), this.options.type instanceof x) {
      if (i.pushCode(`${r} = {};`), i.useContextVariables) {
        const s = i.generateVariable();
        i.pushCode(`${r}.$parent = ${s};`), i.pushCode(`${r}.$root = ${s}.$root;`);
      }
      i.pushPath(this.varName), this.options.type.generate(i), i.popPath(this.varName), i.useContextVariables && (i.pushCode(`delete ${r}.$parent;`), i.pushCode(`delete ${r}.$root;`));
    } else if (Z.has(this.options.type)) {
      const s = i.generateTmpVariable();
      if (i.pushCode(`var ${s} = ${ie + this.options.type}(offset, {`), i.useContextVariables) {
        const o = i.generateVariable();
        i.pushCode(`$parent: ${o},`), i.pushCode(`$root: ${o}.$root,`);
      }
      i.pushCode("});"), i.pushCode(`${r} = ${s}.result; offset = ${s}.offset;`), this.options.type !== this.alias && i.addReference(this.options.type);
    } else if (Object.keys(L).indexOf(this.options.type) >= 0) {
      const s = Ve[e], o = Fe[e];
      i.pushCode(`${r} = dataView.get${s}(offset, ${o});`), i.pushCode(`offset += ${L[e]};`);
    }
    i.pushCode(`offset = ${a};`);
  }
  generateSaveOffset(i) {
    const e = i.generateVariable(this.varName);
    i.pushCode(`${e} = offset`);
  }
}
class X {
  constructor(i, e) {
    this.ranges = arguments.length === 2 ? [{ min: i, max: e }] : 0 in i ? Object.assign({}, i) : [i];
  }
  min() {
    return this.ranges[0].min;
  }
  max() {
    return this.ranges[this.ranges.length - 1].max;
  }
  contains(i) {
    for (let e = 0; e < this.ranges.length; e += 1) {
      const n = this.ranges[e];
      if (n.min <= i && n.max >= i)
        return !0;
    }
    return !1;
  }
  isContiguous() {
    return this.ranges.length > 1;
  }
  getRanges() {
    return this.ranges.map((i) => new X(i.min, i.max));
  }
  toString() {
    return this.ranges.map((i) => `[${i.min}-${i.max}]`).join(",");
  }
  union(i) {
    const e = this.getRanges().concat(i.getRanges()).sort(this.rangeOrder), n = [];
    let a = e[0];
    for (let r = 1; r < e.length; r += 1) {
      const s = e[r];
      s.min() > a.max() + 1 ? (n.push(a), a = s) : s.max() > a.max() && (a = new X(a.min(), s.max()));
    }
    return n.push(a), n.length === 1 ? n[0] : new X(n);
  }
  intersection(i) {
    let e = this, n = i;
    const a = this.ranges(), r = n.ranges(), s = a.length, o = r.length;
    let u = 0, f = 0;
    const l = [];
    for (; u < s && f < o; ) {
      e = a[u], n = r[f];
      const b = Math.max(e.min(), n.min()), d = Math.min(e.max(), n.max());
      d >= b && l.push(new X(b, d)), e.max() > n.max() ? f += 1 : u += 1;
    }
    if (l.length === 0)
      throw new Error("found range of length 0");
    return l.length === 1 ? l[0] : new X(l);
  }
  coverage() {
    let i = 0;
    const e = this.ranges();
    for (let n = 0; n < e.length; n += 1) {
      const a = e[n];
      i += a.max() - a.min() + 1;
    }
    return i;
  }
  rangeOrder(i, e) {
    let n = i, a = e;
    return arguments.length < 2 && (a = n, n = this), n.min() < a.min() ? -1 : n.min() > a.min() ? 1 : n.max() < a.max() ? -1 : a.max() > n.max() ? 1 : 0;
  }
}
/*! pako 2.1.0 https://github.com/nodeca/pako @license (MIT AND Zlib) */
const pn = 4, Dt = 0, Bt = 1, bn = 2;
function pe(t) {
  let i = t.length;
  for (; --i >= 0; )
    t[i] = 0;
}
const gn = 0, $i = 1, wn = 2, mn = 3, yn = 258, Nt = 29, ze = 256, Ce = ze + 1 + Nt, ue = 30, St = 19, Ci = 2 * Ce + 1, ne = 15, tt = 16, vn = 7, Tt = 256, Ei = 16, Ni = 17, Si = 18, gt = (
  /* extra bits for each length code */
  new Uint8Array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0])
), Ge = (
  /* extra bits for each distance code */
  new Uint8Array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13])
), kn = (
  /* extra bits for each bit length code */
  new Uint8Array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7])
), Ti = new Uint8Array([16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]), $n = 512, K = new Array((Ce + 2) * 2);
pe(K);
const ve = new Array(ue * 2);
pe(ve);
const Ee = new Array($n);
pe(Ee);
const Ne = new Array(yn - mn + 1);
pe(Ne);
const xt = new Array(Nt);
pe(xt);
const Ye = new Array(ue);
pe(Ye);
function it(t, i, e, n, a) {
  this.static_tree = t, this.extra_bits = i, this.extra_base = e, this.elems = n, this.max_length = a, this.has_stree = t && t.length;
}
let xi, Ii, Ai;
function nt(t, i) {
  this.dyn_tree = t, this.max_code = 0, this.stat_desc = i;
}
const zi = (t) => t < 256 ? Ee[t] : Ee[256 + (t >>> 7)], Se = (t, i) => {
  t.pending_buf[t.pending++] = i & 255, t.pending_buf[t.pending++] = i >>> 8 & 255;
}, D = (t, i, e) => {
  t.bi_valid > tt - e ? (t.bi_buf |= i << t.bi_valid & 65535, Se(t, t.bi_buf), t.bi_buf = i >> tt - t.bi_valid, t.bi_valid += e - tt) : (t.bi_buf |= i << t.bi_valid & 65535, t.bi_valid += e);
}, M = (t, i, e) => {
  D(
    t,
    e[i * 2],
    e[i * 2 + 1]
    /*.Len*/
  );
}, Oi = (t, i) => {
  let e = 0;
  do
    e |= t & 1, t >>>= 1, e <<= 1;
  while (--i > 0);
  return e >>> 1;
}, Cn = (t) => {
  t.bi_valid === 16 ? (Se(t, t.bi_buf), t.bi_buf = 0, t.bi_valid = 0) : t.bi_valid >= 8 && (t.pending_buf[t.pending++] = t.bi_buf & 255, t.bi_buf >>= 8, t.bi_valid -= 8);
}, En = (t, i) => {
  const e = i.dyn_tree, n = i.max_code, a = i.stat_desc.static_tree, r = i.stat_desc.has_stree, s = i.stat_desc.extra_bits, o = i.stat_desc.extra_base, u = i.stat_desc.max_length;
  let f, l, b, d, h, p, k = 0;
  for (d = 0; d <= ne; d++)
    t.bl_count[d] = 0;
  for (e[t.heap[t.heap_max] * 2 + 1] = 0, f = t.heap_max + 1; f < Ci; f++)
    l = t.heap[f], d = e[e[l * 2 + 1] * 2 + 1] + 1, d > u && (d = u, k++), e[l * 2 + 1] = d, !(l > n) && (t.bl_count[d]++, h = 0, l >= o && (h = s[l - o]), p = e[l * 2], t.opt_len += p * (d + h), r && (t.static_len += p * (a[l * 2 + 1] + h)));
  if (k !== 0) {
    do {
      for (d = u - 1; t.bl_count[d] === 0; )
        d--;
      t.bl_count[d]--, t.bl_count[d + 1] += 2, t.bl_count[u]--, k -= 2;
    } while (k > 0);
    for (d = u; d !== 0; d--)
      for (l = t.bl_count[d]; l !== 0; )
        b = t.heap[--f], !(b > n) && (e[b * 2 + 1] !== d && (t.opt_len += (d - e[b * 2 + 1]) * e[b * 2], e[b * 2 + 1] = d), l--);
  }
}, Ri = (t, i, e) => {
  const n = new Array(ne + 1);
  let a = 0, r, s;
  for (r = 1; r <= ne; r++)
    a = a + e[r - 1] << 1, n[r] = a;
  for (s = 0; s <= i; s++) {
    let o = t[s * 2 + 1];
    o !== 0 && (t[s * 2] = Oi(n[o]++, o));
  }
}, Nn = () => {
  let t, i, e, n, a;
  const r = new Array(ne + 1);
  for (e = 0, n = 0; n < Nt - 1; n++)
    for (xt[n] = e, t = 0; t < 1 << gt[n]; t++)
      Ne[e++] = n;
  for (Ne[e - 1] = n, a = 0, n = 0; n < 16; n++)
    for (Ye[n] = a, t = 0; t < 1 << Ge[n]; t++)
      Ee[a++] = n;
  for (a >>= 7; n < ue; n++)
    for (Ye[n] = a << 7, t = 0; t < 1 << Ge[n] - 7; t++)
      Ee[256 + a++] = n;
  for (i = 0; i <= ne; i++)
    r[i] = 0;
  for (t = 0; t <= 143; )
    K[t * 2 + 1] = 8, t++, r[8]++;
  for (; t <= 255; )
    K[t * 2 + 1] = 9, t++, r[9]++;
  for (; t <= 279; )
    K[t * 2 + 1] = 7, t++, r[7]++;
  for (; t <= 287; )
    K[t * 2 + 1] = 8, t++, r[8]++;
  for (Ri(K, Ce + 1, r), t = 0; t < ue; t++)
    ve[t * 2 + 1] = 5, ve[t * 2] = Oi(t, 5);
  xi = new it(K, gt, ze + 1, Ce, ne), Ii = new it(ve, Ge, 0, ue, ne), Ai = new it(new Array(0), kn, 0, St, vn);
}, Di = (t) => {
  let i;
  for (i = 0; i < Ce; i++)
    t.dyn_ltree[i * 2] = 0;
  for (i = 0; i < ue; i++)
    t.dyn_dtree[i * 2] = 0;
  for (i = 0; i < St; i++)
    t.bl_tree[i * 2] = 0;
  t.dyn_ltree[Tt * 2] = 1, t.opt_len = t.static_len = 0, t.sym_next = t.matches = 0;
}, Bi = (t) => {
  t.bi_valid > 8 ? Se(t, t.bi_buf) : t.bi_valid > 0 && (t.pending_buf[t.pending++] = t.bi_buf), t.bi_buf = 0, t.bi_valid = 0;
}, Ut = (t, i, e, n) => {
  const a = i * 2, r = e * 2;
  return t[a] < t[r] || t[a] === t[r] && n[i] <= n[e];
}, at = (t, i, e) => {
  const n = t.heap[e];
  let a = e << 1;
  for (; a <= t.heap_len && (a < t.heap_len && Ut(i, t.heap[a + 1], t.heap[a], t.depth) && a++, !Ut(i, n, t.heap[a], t.depth)); )
    t.heap[e] = t.heap[a], e = a, a <<= 1;
  t.heap[e] = n;
}, Vt = (t, i, e) => {
  let n, a, r = 0, s, o;
  if (t.sym_next !== 0)
    do
      n = t.pending_buf[t.sym_buf + r++] & 255, n += (t.pending_buf[t.sym_buf + r++] & 255) << 8, a = t.pending_buf[t.sym_buf + r++], n === 0 ? M(t, a, i) : (s = Ne[a], M(t, s + ze + 1, i), o = gt[s], o !== 0 && (a -= xt[s], D(t, a, o)), n--, s = zi(n), M(t, s, e), o = Ge[s], o !== 0 && (n -= Ye[s], D(t, n, o)));
    while (r < t.sym_next);
  M(t, Tt, i);
}, wt = (t, i) => {
  const e = i.dyn_tree, n = i.stat_desc.static_tree, a = i.stat_desc.has_stree, r = i.stat_desc.elems;
  let s, o, u = -1, f;
  for (t.heap_len = 0, t.heap_max = Ci, s = 0; s < r; s++)
    e[s * 2] !== 0 ? (t.heap[++t.heap_len] = u = s, t.depth[s] = 0) : e[s * 2 + 1] = 0;
  for (; t.heap_len < 2; )
    f = t.heap[++t.heap_len] = u < 2 ? ++u : 0, e[f * 2] = 1, t.depth[f] = 0, t.opt_len--, a && (t.static_len -= n[f * 2 + 1]);
  for (i.max_code = u, s = t.heap_len >> 1; s >= 1; s--)
    at(t, e, s);
  f = r;
  do
    s = t.heap[
      1
      /*SMALLEST*/
    ], t.heap[
      1
      /*SMALLEST*/
    ] = t.heap[t.heap_len--], at(
      t,
      e,
      1
      /*SMALLEST*/
    ), o = t.heap[
      1
      /*SMALLEST*/
    ], t.heap[--t.heap_max] = s, t.heap[--t.heap_max] = o, e[f * 2] = e[s * 2] + e[o * 2], t.depth[f] = (t.depth[s] >= t.depth[o] ? t.depth[s] : t.depth[o]) + 1, e[s * 2 + 1] = e[o * 2 + 1] = f, t.heap[
      1
      /*SMALLEST*/
    ] = f++, at(
      t,
      e,
      1
      /*SMALLEST*/
    );
  while (t.heap_len >= 2);
  t.heap[--t.heap_max] = t.heap[
    1
    /*SMALLEST*/
  ], En(t, i), Ri(e, u, t.bl_count);
}, Ft = (t, i, e) => {
  let n, a = -1, r, s = i[0 * 2 + 1], o = 0, u = 7, f = 4;
  for (s === 0 && (u = 138, f = 3), i[(e + 1) * 2 + 1] = 65535, n = 0; n <= e; n++)
    r = s, s = i[(n + 1) * 2 + 1], !(++o < u && r === s) && (o < f ? t.bl_tree[r * 2] += o : r !== 0 ? (r !== a && t.bl_tree[r * 2]++, t.bl_tree[Ei * 2]++) : o <= 10 ? t.bl_tree[Ni * 2]++ : t.bl_tree[Si * 2]++, o = 0, a = r, s === 0 ? (u = 138, f = 3) : r === s ? (u = 6, f = 3) : (u = 7, f = 4));
}, Zt = (t, i, e) => {
  let n, a = -1, r, s = i[0 * 2 + 1], o = 0, u = 7, f = 4;
  for (s === 0 && (u = 138, f = 3), n = 0; n <= e; n++)
    if (r = s, s = i[(n + 1) * 2 + 1], !(++o < u && r === s)) {
      if (o < f)
        do
          M(t, r, t.bl_tree);
        while (--o !== 0);
      else
        r !== 0 ? (r !== a && (M(t, r, t.bl_tree), o--), M(t, Ei, t.bl_tree), D(t, o - 3, 2)) : o <= 10 ? (M(t, Ni, t.bl_tree), D(t, o - 3, 3)) : (M(t, Si, t.bl_tree), D(t, o - 11, 7));
      o = 0, a = r, s === 0 ? (u = 138, f = 3) : r === s ? (u = 6, f = 3) : (u = 7, f = 4);
    }
}, Sn = (t) => {
  let i;
  for (Ft(t, t.dyn_ltree, t.l_desc.max_code), Ft(t, t.dyn_dtree, t.d_desc.max_code), wt(t, t.bl_desc), i = St - 1; i >= 3 && t.bl_tree[Ti[i] * 2 + 1] === 0; i--)
    ;
  return t.opt_len += 3 * (i + 1) + 5 + 5 + 4, i;
}, Tn = (t, i, e, n) => {
  let a;
  for (D(t, i - 257, 5), D(t, e - 1, 5), D(t, n - 4, 4), a = 0; a < n; a++)
    D(t, t.bl_tree[Ti[a] * 2 + 1], 3);
  Zt(t, t.dyn_ltree, i - 1), Zt(t, t.dyn_dtree, e - 1);
}, xn = (t) => {
  let i = 4093624447, e;
  for (e = 0; e <= 31; e++, i >>>= 1)
    if (i & 1 && t.dyn_ltree[e * 2] !== 0)
      return Dt;
  if (t.dyn_ltree[9 * 2] !== 0 || t.dyn_ltree[10 * 2] !== 0 || t.dyn_ltree[13 * 2] !== 0)
    return Bt;
  for (e = 32; e < ze; e++)
    if (t.dyn_ltree[e * 2] !== 0)
      return Bt;
  return Dt;
};
let Lt = !1;
const In = (t) => {
  Lt || (Nn(), Lt = !0), t.l_desc = new nt(t.dyn_ltree, xi), t.d_desc = new nt(t.dyn_dtree, Ii), t.bl_desc = new nt(t.bl_tree, Ai), t.bi_buf = 0, t.bi_valid = 0, Di(t);
}, Ui = (t, i, e, n) => {
  D(t, (gn << 1) + (n ? 1 : 0), 3), Bi(t), Se(t, e), Se(t, ~e), e && t.pending_buf.set(t.window.subarray(i, i + e), t.pending), t.pending += e;
}, An = (t) => {
  D(t, $i << 1, 3), M(t, Tt, K), Cn(t);
}, zn = (t, i, e, n) => {
  let a, r, s = 0;
  t.level > 0 ? (t.strm.data_type === bn && (t.strm.data_type = xn(t)), wt(t, t.l_desc), wt(t, t.d_desc), s = Sn(t), a = t.opt_len + 3 + 7 >>> 3, r = t.static_len + 3 + 7 >>> 3, r <= a && (a = r)) : a = r = e + 5, e + 4 <= a && i !== -1 ? Ui(t, i, e, n) : t.strategy === pn || r === a ? (D(t, ($i << 1) + (n ? 1 : 0), 3), Vt(t, K, ve)) : (D(t, (wn << 1) + (n ? 1 : 0), 3), Tn(t, t.l_desc.max_code + 1, t.d_desc.max_code + 1, s + 1), Vt(t, t.dyn_ltree, t.dyn_dtree)), Di(t), n && Bi(t);
}, On = (t, i, e) => (t.pending_buf[t.sym_buf + t.sym_next++] = i, t.pending_buf[t.sym_buf + t.sym_next++] = i >> 8, t.pending_buf[t.sym_buf + t.sym_next++] = e, i === 0 ? t.dyn_ltree[e * 2]++ : (t.matches++, i--, t.dyn_ltree[(Ne[e] + ze + 1) * 2]++, t.dyn_dtree[zi(i) * 2]++), t.sym_next === t.sym_end);
var Rn = In, Dn = Ui, Bn = zn, Un = On, Vn = An, Fn = {
  _tr_init: Rn,
  _tr_stored_block: Dn,
  _tr_flush_block: Bn,
  _tr_tally: Un,
  _tr_align: Vn
};
const Zn = (t, i, e, n) => {
  let a = t & 65535 | 0, r = t >>> 16 & 65535 | 0, s = 0;
  for (; e !== 0; ) {
    s = e > 2e3 ? 2e3 : e, e -= s;
    do
      a = a + i[n++] | 0, r = r + a | 0;
    while (--s);
    a %= 65521, r %= 65521;
  }
  return a | r << 16 | 0;
};
var Te = Zn;
const Ln = () => {
  let t, i = [];
  for (var e = 0; e < 256; e++) {
    t = e;
    for (var n = 0; n < 8; n++)
      t = t & 1 ? 3988292384 ^ t >>> 1 : t >>> 1;
    i[e] = t;
  }
  return i;
}, Pn = new Uint32Array(Ln()), Mn = (t, i, e, n) => {
  const a = Pn, r = n + e;
  t ^= -1;
  for (let s = n; s < r; s++)
    t = t >>> 8 ^ a[(t ^ i[s]) & 255];
  return t ^ -1;
};
var z = Mn, de = {
  2: "need dictionary",
  /* Z_NEED_DICT       2  */
  1: "stream end",
  /* Z_STREAM_END      1  */
  0: "",
  /* Z_OK              0  */
  "-1": "file error",
  /* Z_ERRNO         (-1) */
  "-2": "stream error",
  /* Z_STREAM_ERROR  (-2) */
  "-3": "data error",
  /* Z_DATA_ERROR    (-3) */
  "-4": "insufficient memory",
  /* Z_MEM_ERROR     (-4) */
  "-5": "buffer error",
  /* Z_BUF_ERROR     (-5) */
  "-6": "incompatible version"
  /* Z_VERSION_ERROR (-6) */
}, Oe = {
  /* Allowed flush values; see deflate() and inflate() below for details */
  Z_NO_FLUSH: 0,
  Z_PARTIAL_FLUSH: 1,
  Z_SYNC_FLUSH: 2,
  Z_FULL_FLUSH: 3,
  Z_FINISH: 4,
  Z_BLOCK: 5,
  Z_TREES: 6,
  /* Return codes for the compression/decompression functions. Negative values
  * are errors, positive values are used for special but normal events.
  */
  Z_OK: 0,
  Z_STREAM_END: 1,
  Z_NEED_DICT: 2,
  Z_ERRNO: -1,
  Z_STREAM_ERROR: -2,
  Z_DATA_ERROR: -3,
  Z_MEM_ERROR: -4,
  Z_BUF_ERROR: -5,
  //Z_VERSION_ERROR: -6,
  /* compression levels */
  Z_NO_COMPRESSION: 0,
  Z_BEST_SPEED: 1,
  Z_BEST_COMPRESSION: 9,
  Z_DEFAULT_COMPRESSION: -1,
  Z_FILTERED: 1,
  Z_HUFFMAN_ONLY: 2,
  Z_RLE: 3,
  Z_FIXED: 4,
  Z_DEFAULT_STRATEGY: 0,
  /* Possible values of the data_type field (though see inflate()) */
  Z_BINARY: 0,
  Z_TEXT: 1,
  //Z_ASCII:                1, // = Z_TEXT (deprecated)
  Z_UNKNOWN: 2,
  /* The deflate compression method */
  Z_DEFLATED: 8
  //Z_NULL:                 null // Use -1 or null inline, depending on var type
};
const { _tr_init: Hn, _tr_stored_block: mt, _tr_flush_block: jn, _tr_tally: J, _tr_align: Gn } = Fn, {
  Z_NO_FLUSH: Q,
  Z_PARTIAL_FLUSH: Kn,
  Z_FULL_FLUSH: Yn,
  Z_FINISH: V,
  Z_BLOCK: Pt,
  Z_OK: O,
  Z_STREAM_END: Mt,
  Z_STREAM_ERROR: H,
  Z_DATA_ERROR: Xn,
  Z_BUF_ERROR: rt,
  Z_DEFAULT_COMPRESSION: Wn,
  Z_FILTERED: qn,
  Z_HUFFMAN_ONLY: Ze,
  Z_RLE: Jn,
  Z_FIXED: Qn,
  Z_DEFAULT_STRATEGY: ea,
  Z_UNKNOWN: ta,
  Z_DEFLATED: Je
} = Oe, ia = 9, na = 15, aa = 8, ra = 29, sa = 256, yt = sa + 1 + ra, oa = 30, fa = 19, la = 2 * yt + 1, ha = 15, T = 3, q = 258, j = q + T + 1, ua = 32, ce = 42, It = 57, vt = 69, kt = 73, $t = 91, Ct = 103, ae = 113, me = 666, R = 1, be = 2, se = 3, ge = 4, da = 3, re = (t, i) => (t.msg = de[i], i), Ht = (t) => t * 2 - (t > 4 ? 9 : 0), W = (t) => {
  let i = t.length;
  for (; --i >= 0; )
    t[i] = 0;
}, ca = (t) => {
  let i, e, n, a = t.w_size;
  i = t.hash_size, n = i;
  do
    e = t.head[--n], t.head[n] = e >= a ? e - a : 0;
  while (--i);
  i = a, n = i;
  do
    e = t.prev[--n], t.prev[n] = e >= a ? e - a : 0;
  while (--i);
};
let _a = (t, i, e) => (i << t.hash_shift ^ e) & t.hash_mask, ee = _a;
const B = (t) => {
  const i = t.state;
  let e = i.pending;
  e > t.avail_out && (e = t.avail_out), e !== 0 && (t.output.set(i.pending_buf.subarray(i.pending_out, i.pending_out + e), t.next_out), t.next_out += e, i.pending_out += e, t.total_out += e, t.avail_out -= e, i.pending -= e, i.pending === 0 && (i.pending_out = 0));
}, U = (t, i) => {
  jn(t, t.block_start >= 0 ? t.block_start : -1, t.strstart - t.block_start, i), t.block_start = t.strstart, B(t.strm);
}, I = (t, i) => {
  t.pending_buf[t.pending++] = i;
}, we = (t, i) => {
  t.pending_buf[t.pending++] = i >>> 8 & 255, t.pending_buf[t.pending++] = i & 255;
}, Et = (t, i, e, n) => {
  let a = t.avail_in;
  return a > n && (a = n), a === 0 ? 0 : (t.avail_in -= a, i.set(t.input.subarray(t.next_in, t.next_in + a), e), t.state.wrap === 1 ? t.adler = Te(t.adler, i, a, e) : t.state.wrap === 2 && (t.adler = z(t.adler, i, a, e)), t.next_in += a, t.total_in += a, a);
}, Vi = (t, i) => {
  let e = t.max_chain_length, n = t.strstart, a, r, s = t.prev_length, o = t.nice_match;
  const u = t.strstart > t.w_size - j ? t.strstart - (t.w_size - j) : 0, f = t.window, l = t.w_mask, b = t.prev, d = t.strstart + q;
  let h = f[n + s - 1], p = f[n + s];
  t.prev_length >= t.good_match && (e >>= 2), o > t.lookahead && (o = t.lookahead);
  do
    if (a = i, !(f[a + s] !== p || f[a + s - 1] !== h || f[a] !== f[n] || f[++a] !== f[n + 1])) {
      n += 2, a++;
      do
        ;
      while (f[++n] === f[++a] && f[++n] === f[++a] && f[++n] === f[++a] && f[++n] === f[++a] && f[++n] === f[++a] && f[++n] === f[++a] && f[++n] === f[++a] && f[++n] === f[++a] && n < d);
      if (r = q - (d - n), n = d - q, r > s) {
        if (t.match_start = i, s = r, r >= o)
          break;
        h = f[n + s - 1], p = f[n + s];
      }
    }
  while ((i = b[i & l]) > u && --e !== 0);
  return s <= t.lookahead ? s : t.lookahead;
}, _e = (t) => {
  const i = t.w_size;
  let e, n, a;
  do {
    if (n = t.window_size - t.lookahead - t.strstart, t.strstart >= i + (i - j) && (t.window.set(t.window.subarray(i, i + i - n), 0), t.match_start -= i, t.strstart -= i, t.block_start -= i, t.insert > t.strstart && (t.insert = t.strstart), ca(t), n += i), t.strm.avail_in === 0)
      break;
    if (e = Et(t.strm, t.window, t.strstart + t.lookahead, n), t.lookahead += e, t.lookahead + t.insert >= T)
      for (a = t.strstart - t.insert, t.ins_h = t.window[a], t.ins_h = ee(t, t.ins_h, t.window[a + 1]); t.insert && (t.ins_h = ee(t, t.ins_h, t.window[a + T - 1]), t.prev[a & t.w_mask] = t.head[t.ins_h], t.head[t.ins_h] = a, a++, t.insert--, !(t.lookahead + t.insert < T)); )
        ;
  } while (t.lookahead < j && t.strm.avail_in !== 0);
}, Fi = (t, i) => {
  let e = t.pending_buf_size - 5 > t.w_size ? t.w_size : t.pending_buf_size - 5, n, a, r, s = 0, o = t.strm.avail_in;
  do {
    if (n = 65535, r = t.bi_valid + 42 >> 3, t.strm.avail_out < r || (r = t.strm.avail_out - r, a = t.strstart - t.block_start, n > a + t.strm.avail_in && (n = a + t.strm.avail_in), n > r && (n = r), n < e && (n === 0 && i !== V || i === Q || n !== a + t.strm.avail_in)))
      break;
    s = i === V && n === a + t.strm.avail_in ? 1 : 0, mt(t, 0, 0, s), t.pending_buf[t.pending - 4] = n, t.pending_buf[t.pending - 3] = n >> 8, t.pending_buf[t.pending - 2] = ~n, t.pending_buf[t.pending - 1] = ~n >> 8, B(t.strm), a && (a > n && (a = n), t.strm.output.set(t.window.subarray(t.block_start, t.block_start + a), t.strm.next_out), t.strm.next_out += a, t.strm.avail_out -= a, t.strm.total_out += a, t.block_start += a, n -= a), n && (Et(t.strm, t.strm.output, t.strm.next_out, n), t.strm.next_out += n, t.strm.avail_out -= n, t.strm.total_out += n);
  } while (s === 0);
  return o -= t.strm.avail_in, o && (o >= t.w_size ? (t.matches = 2, t.window.set(t.strm.input.subarray(t.strm.next_in - t.w_size, t.strm.next_in), 0), t.strstart = t.w_size, t.insert = t.strstart) : (t.window_size - t.strstart <= o && (t.strstart -= t.w_size, t.window.set(t.window.subarray(t.w_size, t.w_size + t.strstart), 0), t.matches < 2 && t.matches++, t.insert > t.strstart && (t.insert = t.strstart)), t.window.set(t.strm.input.subarray(t.strm.next_in - o, t.strm.next_in), t.strstart), t.strstart += o, t.insert += o > t.w_size - t.insert ? t.w_size - t.insert : o), t.block_start = t.strstart), t.high_water < t.strstart && (t.high_water = t.strstart), s ? ge : i !== Q && i !== V && t.strm.avail_in === 0 && t.strstart === t.block_start ? be : (r = t.window_size - t.strstart, t.strm.avail_in > r && t.block_start >= t.w_size && (t.block_start -= t.w_size, t.strstart -= t.w_size, t.window.set(t.window.subarray(t.w_size, t.w_size + t.strstart), 0), t.matches < 2 && t.matches++, r += t.w_size, t.insert > t.strstart && (t.insert = t.strstart)), r > t.strm.avail_in && (r = t.strm.avail_in), r && (Et(t.strm, t.window, t.strstart, r), t.strstart += r, t.insert += r > t.w_size - t.insert ? t.w_size - t.insert : r), t.high_water < t.strstart && (t.high_water = t.strstart), r = t.bi_valid + 42 >> 3, r = t.pending_buf_size - r > 65535 ? 65535 : t.pending_buf_size - r, e = r > t.w_size ? t.w_size : r, a = t.strstart - t.block_start, (a >= e || (a || i === V) && i !== Q && t.strm.avail_in === 0 && a <= r) && (n = a > r ? r : a, s = i === V && t.strm.avail_in === 0 && n === a ? 1 : 0, mt(t, t.block_start, n, s), t.block_start += n, B(t.strm)), s ? se : R);
}, st = (t, i) => {
  let e, n;
  for (; ; ) {
    if (t.lookahead < j) {
      if (_e(t), t.lookahead < j && i === Q)
        return R;
      if (t.lookahead === 0)
        break;
    }
    if (e = 0, t.lookahead >= T && (t.ins_h = ee(t, t.ins_h, t.window[t.strstart + T - 1]), e = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h], t.head[t.ins_h] = t.strstart), e !== 0 && t.strstart - e <= t.w_size - j && (t.match_length = Vi(t, e)), t.match_length >= T)
      if (n = J(t, t.strstart - t.match_start, t.match_length - T), t.lookahead -= t.match_length, t.match_length <= t.max_lazy_match && t.lookahead >= T) {
        t.match_length--;
        do
          t.strstart++, t.ins_h = ee(t, t.ins_h, t.window[t.strstart + T - 1]), e = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h], t.head[t.ins_h] = t.strstart;
        while (--t.match_length !== 0);
        t.strstart++;
      } else
        t.strstart += t.match_length, t.match_length = 0, t.ins_h = t.window[t.strstart], t.ins_h = ee(t, t.ins_h, t.window[t.strstart + 1]);
    else
      n = J(t, 0, t.window[t.strstart]), t.lookahead--, t.strstart++;
    if (n && (U(t, !1), t.strm.avail_out === 0))
      return R;
  }
  return t.insert = t.strstart < T - 1 ? t.strstart : T - 1, i === V ? (U(t, !0), t.strm.avail_out === 0 ? se : ge) : t.sym_next && (U(t, !1), t.strm.avail_out === 0) ? R : be;
}, le = (t, i) => {
  let e, n, a;
  for (; ; ) {
    if (t.lookahead < j) {
      if (_e(t), t.lookahead < j && i === Q)
        return R;
      if (t.lookahead === 0)
        break;
    }
    if (e = 0, t.lookahead >= T && (t.ins_h = ee(t, t.ins_h, t.window[t.strstart + T - 1]), e = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h], t.head[t.ins_h] = t.strstart), t.prev_length = t.match_length, t.prev_match = t.match_start, t.match_length = T - 1, e !== 0 && t.prev_length < t.max_lazy_match && t.strstart - e <= t.w_size - j && (t.match_length = Vi(t, e), t.match_length <= 5 && (t.strategy === qn || t.match_length === T && t.strstart - t.match_start > 4096) && (t.match_length = T - 1)), t.prev_length >= T && t.match_length <= t.prev_length) {
      a = t.strstart + t.lookahead - T, n = J(t, t.strstart - 1 - t.prev_match, t.prev_length - T), t.lookahead -= t.prev_length - 1, t.prev_length -= 2;
      do
        ++t.strstart <= a && (t.ins_h = ee(t, t.ins_h, t.window[t.strstart + T - 1]), e = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h], t.head[t.ins_h] = t.strstart);
      while (--t.prev_length !== 0);
      if (t.match_available = 0, t.match_length = T - 1, t.strstart++, n && (U(t, !1), t.strm.avail_out === 0))
        return R;
    } else if (t.match_available) {
      if (n = J(t, 0, t.window[t.strstart - 1]), n && U(t, !1), t.strstart++, t.lookahead--, t.strm.avail_out === 0)
        return R;
    } else
      t.match_available = 1, t.strstart++, t.lookahead--;
  }
  return t.match_available && (n = J(t, 0, t.window[t.strstart - 1]), t.match_available = 0), t.insert = t.strstart < T - 1 ? t.strstart : T - 1, i === V ? (U(t, !0), t.strm.avail_out === 0 ? se : ge) : t.sym_next && (U(t, !1), t.strm.avail_out === 0) ? R : be;
}, pa = (t, i) => {
  let e, n, a, r;
  const s = t.window;
  for (; ; ) {
    if (t.lookahead <= q) {
      if (_e(t), t.lookahead <= q && i === Q)
        return R;
      if (t.lookahead === 0)
        break;
    }
    if (t.match_length = 0, t.lookahead >= T && t.strstart > 0 && (a = t.strstart - 1, n = s[a], n === s[++a] && n === s[++a] && n === s[++a])) {
      r = t.strstart + q;
      do
        ;
      while (n === s[++a] && n === s[++a] && n === s[++a] && n === s[++a] && n === s[++a] && n === s[++a] && n === s[++a] && n === s[++a] && a < r);
      t.match_length = q - (r - a), t.match_length > t.lookahead && (t.match_length = t.lookahead);
    }
    if (t.match_length >= T ? (e = J(t, 1, t.match_length - T), t.lookahead -= t.match_length, t.strstart += t.match_length, t.match_length = 0) : (e = J(t, 0, t.window[t.strstart]), t.lookahead--, t.strstart++), e && (U(t, !1), t.strm.avail_out === 0))
      return R;
  }
  return t.insert = 0, i === V ? (U(t, !0), t.strm.avail_out === 0 ? se : ge) : t.sym_next && (U(t, !1), t.strm.avail_out === 0) ? R : be;
}, ba = (t, i) => {
  let e;
  for (; ; ) {
    if (t.lookahead === 0 && (_e(t), t.lookahead === 0)) {
      if (i === Q)
        return R;
      break;
    }
    if (t.match_length = 0, e = J(t, 0, t.window[t.strstart]), t.lookahead--, t.strstart++, e && (U(t, !1), t.strm.avail_out === 0))
      return R;
  }
  return t.insert = 0, i === V ? (U(t, !0), t.strm.avail_out === 0 ? se : ge) : t.sym_next && (U(t, !1), t.strm.avail_out === 0) ? R : be;
};
function P(t, i, e, n, a) {
  this.good_length = t, this.max_lazy = i, this.nice_length = e, this.max_chain = n, this.func = a;
}
const ye = [
  /*      good lazy nice chain */
  new P(0, 0, 0, 0, Fi),
  /* 0 store only */
  new P(4, 4, 8, 4, st),
  /* 1 max speed, no lazy matches */
  new P(4, 5, 16, 8, st),
  /* 2 */
  new P(4, 6, 32, 32, st),
  /* 3 */
  new P(4, 4, 16, 16, le),
  /* 4 lazy matches */
  new P(8, 16, 32, 32, le),
  /* 5 */
  new P(8, 16, 128, 128, le),
  /* 6 */
  new P(8, 32, 128, 256, le),
  /* 7 */
  new P(32, 128, 258, 1024, le),
  /* 8 */
  new P(32, 258, 258, 4096, le)
  /* 9 max compression */
], ga = (t) => {
  t.window_size = 2 * t.w_size, W(t.head), t.max_lazy_match = ye[t.level].max_lazy, t.good_match = ye[t.level].good_length, t.nice_match = ye[t.level].nice_length, t.max_chain_length = ye[t.level].max_chain, t.strstart = 0, t.block_start = 0, t.lookahead = 0, t.insert = 0, t.match_length = t.prev_length = T - 1, t.match_available = 0, t.ins_h = 0;
};
function wa() {
  this.strm = null, this.status = 0, this.pending_buf = null, this.pending_buf_size = 0, this.pending_out = 0, this.pending = 0, this.wrap = 0, this.gzhead = null, this.gzindex = 0, this.method = Je, this.last_flush = -1, this.w_size = 0, this.w_bits = 0, this.w_mask = 0, this.window = null, this.window_size = 0, this.prev = null, this.head = null, this.ins_h = 0, this.hash_size = 0, this.hash_bits = 0, this.hash_mask = 0, this.hash_shift = 0, this.block_start = 0, this.match_length = 0, this.prev_match = 0, this.match_available = 0, this.strstart = 0, this.match_start = 0, this.lookahead = 0, this.prev_length = 0, this.max_chain_length = 0, this.max_lazy_match = 0, this.level = 0, this.strategy = 0, this.good_match = 0, this.nice_match = 0, this.dyn_ltree = new Uint16Array(la * 2), this.dyn_dtree = new Uint16Array((2 * oa + 1) * 2), this.bl_tree = new Uint16Array((2 * fa + 1) * 2), W(this.dyn_ltree), W(this.dyn_dtree), W(this.bl_tree), this.l_desc = null, this.d_desc = null, this.bl_desc = null, this.bl_count = new Uint16Array(ha + 1), this.heap = new Uint16Array(2 * yt + 1), W(this.heap), this.heap_len = 0, this.heap_max = 0, this.depth = new Uint16Array(2 * yt + 1), W(this.depth), this.sym_buf = 0, this.lit_bufsize = 0, this.sym_next = 0, this.sym_end = 0, this.opt_len = 0, this.static_len = 0, this.matches = 0, this.insert = 0, this.bi_buf = 0, this.bi_valid = 0;
}
const Re = (t) => {
  if (!t)
    return 1;
  const i = t.state;
  return !i || i.strm !== t || i.status !== ce && //#ifdef GZIP
  i.status !== It && //#endif
  i.status !== vt && i.status !== kt && i.status !== $t && i.status !== Ct && i.status !== ae && i.status !== me ? 1 : 0;
}, Zi = (t) => {
  if (Re(t))
    return re(t, H);
  t.total_in = t.total_out = 0, t.data_type = ta;
  const i = t.state;
  return i.pending = 0, i.pending_out = 0, i.wrap < 0 && (i.wrap = -i.wrap), i.status = //#ifdef GZIP
  i.wrap === 2 ? It : (
    //#endif
    i.wrap ? ce : ae
  ), t.adler = i.wrap === 2 ? 0 : 1, i.last_flush = -2, Hn(i), O;
}, Li = (t) => {
  const i = Zi(t);
  return i === O && ga(t.state), i;
}, ma = (t, i) => Re(t) || t.state.wrap !== 2 ? H : (t.state.gzhead = i, O), Pi = (t, i, e, n, a, r) => {
  if (!t)
    return H;
  let s = 1;
  if (i === Wn && (i = 6), n < 0 ? (s = 0, n = -n) : n > 15 && (s = 2, n -= 16), a < 1 || a > ia || e !== Je || n < 8 || n > 15 || i < 0 || i > 9 || r < 0 || r > Qn || n === 8 && s !== 1)
    return re(t, H);
  n === 8 && (n = 9);
  const o = new wa();
  return t.state = o, o.strm = t, o.status = ce, o.wrap = s, o.gzhead = null, o.w_bits = n, o.w_size = 1 << o.w_bits, o.w_mask = o.w_size - 1, o.hash_bits = a + 7, o.hash_size = 1 << o.hash_bits, o.hash_mask = o.hash_size - 1, o.hash_shift = ~~((o.hash_bits + T - 1) / T), o.window = new Uint8Array(o.w_size * 2), o.head = new Uint16Array(o.hash_size), o.prev = new Uint16Array(o.w_size), o.lit_bufsize = 1 << a + 6, o.pending_buf_size = o.lit_bufsize * 4, o.pending_buf = new Uint8Array(o.pending_buf_size), o.sym_buf = o.lit_bufsize, o.sym_end = (o.lit_bufsize - 1) * 3, o.level = i, o.strategy = r, o.method = e, Li(t);
}, ya = (t, i) => Pi(t, i, Je, na, aa, ea), va = (t, i) => {
  if (Re(t) || i > Pt || i < 0)
    return t ? re(t, H) : H;
  const e = t.state;
  if (!t.output || t.avail_in !== 0 && !t.input || e.status === me && i !== V)
    return re(t, t.avail_out === 0 ? rt : H);
  const n = e.last_flush;
  if (e.last_flush = i, e.pending !== 0) {
    if (B(t), t.avail_out === 0)
      return e.last_flush = -1, O;
  } else if (t.avail_in === 0 && Ht(i) <= Ht(n) && i !== V)
    return re(t, rt);
  if (e.status === me && t.avail_in !== 0)
    return re(t, rt);
  if (e.status === ce && e.wrap === 0 && (e.status = ae), e.status === ce) {
    let a = Je + (e.w_bits - 8 << 4) << 8, r = -1;
    if (e.strategy >= Ze || e.level < 2 ? r = 0 : e.level < 6 ? r = 1 : e.level === 6 ? r = 2 : r = 3, a |= r << 6, e.strstart !== 0 && (a |= ua), a += 31 - a % 31, we(e, a), e.strstart !== 0 && (we(e, t.adler >>> 16), we(e, t.adler & 65535)), t.adler = 1, e.status = ae, B(t), e.pending !== 0)
      return e.last_flush = -1, O;
  }
  if (e.status === It) {
    if (t.adler = 0, I(e, 31), I(e, 139), I(e, 8), e.gzhead)
      I(
        e,
        (e.gzhead.text ? 1 : 0) + (e.gzhead.hcrc ? 2 : 0) + (e.gzhead.extra ? 4 : 0) + (e.gzhead.name ? 8 : 0) + (e.gzhead.comment ? 16 : 0)
      ), I(e, e.gzhead.time & 255), I(e, e.gzhead.time >> 8 & 255), I(e, e.gzhead.time >> 16 & 255), I(e, e.gzhead.time >> 24 & 255), I(e, e.level === 9 ? 2 : e.strategy >= Ze || e.level < 2 ? 4 : 0), I(e, e.gzhead.os & 255), e.gzhead.extra && e.gzhead.extra.length && (I(e, e.gzhead.extra.length & 255), I(e, e.gzhead.extra.length >> 8 & 255)), e.gzhead.hcrc && (t.adler = z(t.adler, e.pending_buf, e.pending, 0)), e.gzindex = 0, e.status = vt;
    else if (I(e, 0), I(e, 0), I(e, 0), I(e, 0), I(e, 0), I(e, e.level === 9 ? 2 : e.strategy >= Ze || e.level < 2 ? 4 : 0), I(e, da), e.status = ae, B(t), e.pending !== 0)
      return e.last_flush = -1, O;
  }
  if (e.status === vt) {
    if (e.gzhead.extra) {
      let a = e.pending, r = (e.gzhead.extra.length & 65535) - e.gzindex;
      for (; e.pending + r > e.pending_buf_size; ) {
        let o = e.pending_buf_size - e.pending;
        if (e.pending_buf.set(e.gzhead.extra.subarray(e.gzindex, e.gzindex + o), e.pending), e.pending = e.pending_buf_size, e.gzhead.hcrc && e.pending > a && (t.adler = z(t.adler, e.pending_buf, e.pending - a, a)), e.gzindex += o, B(t), e.pending !== 0)
          return e.last_flush = -1, O;
        a = 0, r -= o;
      }
      let s = new Uint8Array(e.gzhead.extra);
      e.pending_buf.set(s.subarray(e.gzindex, e.gzindex + r), e.pending), e.pending += r, e.gzhead.hcrc && e.pending > a && (t.adler = z(t.adler, e.pending_buf, e.pending - a, a)), e.gzindex = 0;
    }
    e.status = kt;
  }
  if (e.status === kt) {
    if (e.gzhead.name) {
      let a = e.pending, r;
      do {
        if (e.pending === e.pending_buf_size) {
          if (e.gzhead.hcrc && e.pending > a && (t.adler = z(t.adler, e.pending_buf, e.pending - a, a)), B(t), e.pending !== 0)
            return e.last_flush = -1, O;
          a = 0;
        }
        e.gzindex < e.gzhead.name.length ? r = e.gzhead.name.charCodeAt(e.gzindex++) & 255 : r = 0, I(e, r);
      } while (r !== 0);
      e.gzhead.hcrc && e.pending > a && (t.adler = z(t.adler, e.pending_buf, e.pending - a, a)), e.gzindex = 0;
    }
    e.status = $t;
  }
  if (e.status === $t) {
    if (e.gzhead.comment) {
      let a = e.pending, r;
      do {
        if (e.pending === e.pending_buf_size) {
          if (e.gzhead.hcrc && e.pending > a && (t.adler = z(t.adler, e.pending_buf, e.pending - a, a)), B(t), e.pending !== 0)
            return e.last_flush = -1, O;
          a = 0;
        }
        e.gzindex < e.gzhead.comment.length ? r = e.gzhead.comment.charCodeAt(e.gzindex++) & 255 : r = 0, I(e, r);
      } while (r !== 0);
      e.gzhead.hcrc && e.pending > a && (t.adler = z(t.adler, e.pending_buf, e.pending - a, a));
    }
    e.status = Ct;
  }
  if (e.status === Ct) {
    if (e.gzhead.hcrc) {
      if (e.pending + 2 > e.pending_buf_size && (B(t), e.pending !== 0))
        return e.last_flush = -1, O;
      I(e, t.adler & 255), I(e, t.adler >> 8 & 255), t.adler = 0;
    }
    if (e.status = ae, B(t), e.pending !== 0)
      return e.last_flush = -1, O;
  }
  if (t.avail_in !== 0 || e.lookahead !== 0 || i !== Q && e.status !== me) {
    let a = e.level === 0 ? Fi(e, i) : e.strategy === Ze ? ba(e, i) : e.strategy === Jn ? pa(e, i) : ye[e.level].func(e, i);
    if ((a === se || a === ge) && (e.status = me), a === R || a === se)
      return t.avail_out === 0 && (e.last_flush = -1), O;
    if (a === be && (i === Kn ? Gn(e) : i !== Pt && (mt(e, 0, 0, !1), i === Yn && (W(e.head), e.lookahead === 0 && (e.strstart = 0, e.block_start = 0, e.insert = 0))), B(t), t.avail_out === 0))
      return e.last_flush = -1, O;
  }
  return i !== V ? O : e.wrap <= 0 ? Mt : (e.wrap === 2 ? (I(e, t.adler & 255), I(e, t.adler >> 8 & 255), I(e, t.adler >> 16 & 255), I(e, t.adler >> 24 & 255), I(e, t.total_in & 255), I(e, t.total_in >> 8 & 255), I(e, t.total_in >> 16 & 255), I(e, t.total_in >> 24 & 255)) : (we(e, t.adler >>> 16), we(e, t.adler & 65535)), B(t), e.wrap > 0 && (e.wrap = -e.wrap), e.pending !== 0 ? O : Mt);
}, ka = (t) => {
  if (Re(t))
    return H;
  const i = t.state.status;
  return t.state = null, i === ae ? re(t, Xn) : O;
}, $a = (t, i) => {
  let e = i.length;
  if (Re(t))
    return H;
  const n = t.state, a = n.wrap;
  if (a === 2 || a === 1 && n.status !== ce || n.lookahead)
    return H;
  if (a === 1 && (t.adler = Te(t.adler, i, e, 0)), n.wrap = 0, e >= n.w_size) {
    a === 0 && (W(n.head), n.strstart = 0, n.block_start = 0, n.insert = 0);
    let u = new Uint8Array(n.w_size);
    u.set(i.subarray(e - n.w_size, e), 0), i = u, e = n.w_size;
  }
  const r = t.avail_in, s = t.next_in, o = t.input;
  for (t.avail_in = e, t.next_in = 0, t.input = i, _e(n); n.lookahead >= T; ) {
    let u = n.strstart, f = n.lookahead - (T - 1);
    do
      n.ins_h = ee(n, n.ins_h, n.window[u + T - 1]), n.prev[u & n.w_mask] = n.head[n.ins_h], n.head[n.ins_h] = u, u++;
    while (--f);
    n.strstart = u, n.lookahead = T - 1, _e(n);
  }
  return n.strstart += n.lookahead, n.block_start = n.strstart, n.insert = n.lookahead, n.lookahead = 0, n.match_length = n.prev_length = T - 1, n.match_available = 0, t.next_in = s, t.input = o, t.avail_in = r, n.wrap = a, O;
};
var Ca = ya, Ea = Pi, Na = Li, Sa = Zi, Ta = ma, xa = va, Ia = ka, Aa = $a, za = "pako deflate (from Nodeca project)", ke = {
  deflateInit: Ca,
  deflateInit2: Ea,
  deflateReset: Na,
  deflateResetKeep: Sa,
  deflateSetHeader: Ta,
  deflate: xa,
  deflateEnd: Ia,
  deflateSetDictionary: Aa,
  deflateInfo: za
};
const Oa = (t, i) => Object.prototype.hasOwnProperty.call(t, i);
var Ra = function(t) {
  const i = Array.prototype.slice.call(arguments, 1);
  for (; i.length; ) {
    const e = i.shift();
    if (e) {
      if (typeof e != "object")
        throw new TypeError(e + "must be non-object");
      for (const n in e)
        Oa(e, n) && (t[n] = e[n]);
    }
  }
  return t;
}, Da = (t) => {
  let i = 0;
  for (let n = 0, a = t.length; n < a; n++)
    i += t[n].length;
  const e = new Uint8Array(i);
  for (let n = 0, a = 0, r = t.length; n < r; n++) {
    let s = t[n];
    e.set(s, a), a += s.length;
  }
  return e;
}, Qe = {
  assign: Ra,
  flattenChunks: Da
};
let Mi = !0;
try {
  String.fromCharCode.apply(null, new Uint8Array(1));
} catch {
  Mi = !1;
}
const xe = new Uint8Array(256);
for (let t = 0; t < 256; t++)
  xe[t] = t >= 252 ? 6 : t >= 248 ? 5 : t >= 240 ? 4 : t >= 224 ? 3 : t >= 192 ? 2 : 1;
xe[254] = xe[254] = 1;
var Ba = (t) => {
  if (typeof TextEncoder == "function" && TextEncoder.prototype.encode)
    return new TextEncoder().encode(t);
  let i, e, n, a, r, s = t.length, o = 0;
  for (a = 0; a < s; a++)
    e = t.charCodeAt(a), (e & 64512) === 55296 && a + 1 < s && (n = t.charCodeAt(a + 1), (n & 64512) === 56320 && (e = 65536 + (e - 55296 << 10) + (n - 56320), a++)), o += e < 128 ? 1 : e < 2048 ? 2 : e < 65536 ? 3 : 4;
  for (i = new Uint8Array(o), r = 0, a = 0; r < o; a++)
    e = t.charCodeAt(a), (e & 64512) === 55296 && a + 1 < s && (n = t.charCodeAt(a + 1), (n & 64512) === 56320 && (e = 65536 + (e - 55296 << 10) + (n - 56320), a++)), e < 128 ? i[r++] = e : e < 2048 ? (i[r++] = 192 | e >>> 6, i[r++] = 128 | e & 63) : e < 65536 ? (i[r++] = 224 | e >>> 12, i[r++] = 128 | e >>> 6 & 63, i[r++] = 128 | e & 63) : (i[r++] = 240 | e >>> 18, i[r++] = 128 | e >>> 12 & 63, i[r++] = 128 | e >>> 6 & 63, i[r++] = 128 | e & 63);
  return i;
};
const Ua = (t, i) => {
  if (i < 65534 && t.subarray && Mi)
    return String.fromCharCode.apply(null, t.length === i ? t : t.subarray(0, i));
  let e = "";
  for (let n = 0; n < i; n++)
    e += String.fromCharCode(t[n]);
  return e;
};
var Va = (t, i) => {
  const e = i || t.length;
  if (typeof TextDecoder == "function" && TextDecoder.prototype.decode)
    return new TextDecoder().decode(t.subarray(0, i));
  let n, a;
  const r = new Array(e * 2);
  for (a = 0, n = 0; n < e; ) {
    let s = t[n++];
    if (s < 128) {
      r[a++] = s;
      continue;
    }
    let o = xe[s];
    if (o > 4) {
      r[a++] = 65533, n += o - 1;
      continue;
    }
    for (s &= o === 2 ? 31 : o === 3 ? 15 : 7; o > 1 && n < e; )
      s = s << 6 | t[n++] & 63, o--;
    if (o > 1) {
      r[a++] = 65533;
      continue;
    }
    s < 65536 ? r[a++] = s : (s -= 65536, r[a++] = 55296 | s >> 10 & 1023, r[a++] = 56320 | s & 1023);
  }
  return Ua(r, a);
}, Fa = (t, i) => {
  i = i || t.length, i > t.length && (i = t.length);
  let e = i - 1;
  for (; e >= 0 && (t[e] & 192) === 128; )
    e--;
  return e < 0 || e === 0 ? i : e + xe[t[e]] > i ? e : i;
}, Ie = {
  string2buf: Ba,
  buf2string: Va,
  utf8border: Fa
};
function Za() {
  this.input = null, this.next_in = 0, this.avail_in = 0, this.total_in = 0, this.output = null, this.next_out = 0, this.avail_out = 0, this.total_out = 0, this.msg = "", this.state = null, this.data_type = 2, this.adler = 0;
}
var Hi = Za;
const ji = Object.prototype.toString, {
  Z_NO_FLUSH: La,
  Z_SYNC_FLUSH: Pa,
  Z_FULL_FLUSH: Ma,
  Z_FINISH: Ha,
  Z_OK: Xe,
  Z_STREAM_END: ja,
  Z_DEFAULT_COMPRESSION: Ga,
  Z_DEFAULT_STRATEGY: Ka,
  Z_DEFLATED: Ya
} = Oe;
function At(t) {
  this.options = Qe.assign({
    level: Ga,
    method: Ya,
    chunkSize: 16384,
    windowBits: 15,
    memLevel: 8,
    strategy: Ka
  }, t || {});
  let i = this.options;
  i.raw && i.windowBits > 0 ? i.windowBits = -i.windowBits : i.gzip && i.windowBits > 0 && i.windowBits < 16 && (i.windowBits += 16), this.err = 0, this.msg = "", this.ended = !1, this.chunks = [], this.strm = new Hi(), this.strm.avail_out = 0;
  let e = ke.deflateInit2(
    this.strm,
    i.level,
    i.method,
    i.windowBits,
    i.memLevel,
    i.strategy
  );
  if (e !== Xe)
    throw new Error(de[e]);
  if (i.header && ke.deflateSetHeader(this.strm, i.header), i.dictionary) {
    let n;
    if (typeof i.dictionary == "string" ? n = Ie.string2buf(i.dictionary) : ji.call(i.dictionary) === "[object ArrayBuffer]" ? n = new Uint8Array(i.dictionary) : n = i.dictionary, e = ke.deflateSetDictionary(this.strm, n), e !== Xe)
      throw new Error(de[e]);
    this._dict_set = !0;
  }
}
At.prototype.push = function(t, i) {
  const e = this.strm, n = this.options.chunkSize;
  let a, r;
  if (this.ended)
    return !1;
  for (i === ~~i ? r = i : r = i === !0 ? Ha : La, typeof t == "string" ? e.input = Ie.string2buf(t) : ji.call(t) === "[object ArrayBuffer]" ? e.input = new Uint8Array(t) : e.input = t, e.next_in = 0, e.avail_in = e.input.length; ; ) {
    if (e.avail_out === 0 && (e.output = new Uint8Array(n), e.next_out = 0, e.avail_out = n), (r === Pa || r === Ma) && e.avail_out <= 6) {
      this.onData(e.output.subarray(0, e.next_out)), e.avail_out = 0;
      continue;
    }
    if (a = ke.deflate(e, r), a === ja)
      return e.next_out > 0 && this.onData(e.output.subarray(0, e.next_out)), a = ke.deflateEnd(this.strm), this.onEnd(a), this.ended = !0, a === Xe;
    if (e.avail_out === 0) {
      this.onData(e.output);
      continue;
    }
    if (r > 0 && e.next_out > 0) {
      this.onData(e.output.subarray(0, e.next_out)), e.avail_out = 0;
      continue;
    }
    if (e.avail_in === 0)
      break;
  }
  return !0;
};
At.prototype.onData = function(t) {
  this.chunks.push(t);
};
At.prototype.onEnd = function(t) {
  t === Xe && (this.result = Qe.flattenChunks(this.chunks)), this.chunks = [], this.err = t, this.msg = this.strm.msg;
};
const Le = 16209, Xa = 16191;
var Wa = function(i, e) {
  let n, a, r, s, o, u, f, l, b, d, h, p, k, C, v, E, $, c, y, S, _, w, m, g;
  const N = i.state;
  n = i.next_in, m = i.input, a = n + (i.avail_in - 5), r = i.next_out, g = i.output, s = r - (e - i.avail_out), o = r + (i.avail_out - 257), u = N.dmax, f = N.wsize, l = N.whave, b = N.wnext, d = N.window, h = N.hold, p = N.bits, k = N.lencode, C = N.distcode, v = (1 << N.lenbits) - 1, E = (1 << N.distbits) - 1;
  e:
    do {
      p < 15 && (h += m[n++] << p, p += 8, h += m[n++] << p, p += 8), $ = k[h & v];
      t:
        for (; ; ) {
          if (c = $ >>> 24, h >>>= c, p -= c, c = $ >>> 16 & 255, c === 0)
            g[r++] = $ & 65535;
          else if (c & 16) {
            y = $ & 65535, c &= 15, c && (p < c && (h += m[n++] << p, p += 8), y += h & (1 << c) - 1, h >>>= c, p -= c), p < 15 && (h += m[n++] << p, p += 8, h += m[n++] << p, p += 8), $ = C[h & E];
            i:
              for (; ; ) {
                if (c = $ >>> 24, h >>>= c, p -= c, c = $ >>> 16 & 255, c & 16) {
                  if (S = $ & 65535, c &= 15, p < c && (h += m[n++] << p, p += 8, p < c && (h += m[n++] << p, p += 8)), S += h & (1 << c) - 1, S > u) {
                    i.msg = "invalid distance too far back", N.mode = Le;
                    break e;
                  }
                  if (h >>>= c, p -= c, c = r - s, S > c) {
                    if (c = S - c, c > l && N.sane) {
                      i.msg = "invalid distance too far back", N.mode = Le;
                      break e;
                    }
                    if (_ = 0, w = d, b === 0) {
                      if (_ += f - c, c < y) {
                        y -= c;
                        do
                          g[r++] = d[_++];
                        while (--c);
                        _ = r - S, w = g;
                      }
                    } else if (b < c) {
                      if (_ += f + b - c, c -= b, c < y) {
                        y -= c;
                        do
                          g[r++] = d[_++];
                        while (--c);
                        if (_ = 0, b < y) {
                          c = b, y -= c;
                          do
                            g[r++] = d[_++];
                          while (--c);
                          _ = r - S, w = g;
                        }
                      }
                    } else if (_ += b - c, c < y) {
                      y -= c;
                      do
                        g[r++] = d[_++];
                      while (--c);
                      _ = r - S, w = g;
                    }
                    for (; y > 2; )
                      g[r++] = w[_++], g[r++] = w[_++], g[r++] = w[_++], y -= 3;
                    y && (g[r++] = w[_++], y > 1 && (g[r++] = w[_++]));
                  } else {
                    _ = r - S;
                    do
                      g[r++] = g[_++], g[r++] = g[_++], g[r++] = g[_++], y -= 3;
                    while (y > 2);
                    y && (g[r++] = g[_++], y > 1 && (g[r++] = g[_++]));
                  }
                } else if (c & 64) {
                  i.msg = "invalid distance code", N.mode = Le;
                  break e;
                } else {
                  $ = C[($ & 65535) + (h & (1 << c) - 1)];
                  continue i;
                }
                break;
              }
          } else if (c & 64)
            if (c & 32) {
              N.mode = Xa;
              break e;
            } else {
              i.msg = "invalid literal/length code", N.mode = Le;
              break e;
            }
          else {
            $ = k[($ & 65535) + (h & (1 << c) - 1)];
            continue t;
          }
          break;
        }
    } while (n < a && r < o);
  y = p >> 3, n -= y, p -= y << 3, h &= (1 << p) - 1, i.next_in = n, i.next_out = r, i.avail_in = n < a ? 5 + (a - n) : 5 - (n - a), i.avail_out = r < o ? 257 + (o - r) : 257 - (r - o), N.hold = h, N.bits = p;
};
const he = 15, jt = 852, Gt = 592, Kt = 0, ot = 1, Yt = 2, qa = new Uint16Array([
  /* Length codes 257..285 base */
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  13,
  15,
  17,
  19,
  23,
  27,
  31,
  35,
  43,
  51,
  59,
  67,
  83,
  99,
  115,
  131,
  163,
  195,
  227,
  258,
  0,
  0
]), Ja = new Uint8Array([
  /* Length codes 257..285 extra */
  16,
  16,
  16,
  16,
  16,
  16,
  16,
  16,
  17,
  17,
  17,
  17,
  18,
  18,
  18,
  18,
  19,
  19,
  19,
  19,
  20,
  20,
  20,
  20,
  21,
  21,
  21,
  21,
  16,
  72,
  78
]), Qa = new Uint16Array([
  /* Distance codes 0..29 base */
  1,
  2,
  3,
  4,
  5,
  7,
  9,
  13,
  17,
  25,
  33,
  49,
  65,
  97,
  129,
  193,
  257,
  385,
  513,
  769,
  1025,
  1537,
  2049,
  3073,
  4097,
  6145,
  8193,
  12289,
  16385,
  24577,
  0,
  0
]), er = new Uint8Array([
  /* Distance codes 0..29 extra */
  16,
  16,
  16,
  16,
  17,
  17,
  18,
  18,
  19,
  19,
  20,
  20,
  21,
  21,
  22,
  22,
  23,
  23,
  24,
  24,
  25,
  25,
  26,
  26,
  27,
  27,
  28,
  28,
  29,
  29,
  64,
  64
]), tr = (t, i, e, n, a, r, s, o) => {
  const u = o.bits;
  let f = 0, l = 0, b = 0, d = 0, h = 0, p = 0, k = 0, C = 0, v = 0, E = 0, $, c, y, S, _, w = null, m;
  const g = new Uint16Array(he + 1), N = new Uint16Array(he + 1);
  let te = null, Ot, Be, Ue;
  for (f = 0; f <= he; f++)
    g[f] = 0;
  for (l = 0; l < n; l++)
    g[i[e + l]]++;
  for (h = u, d = he; d >= 1 && g[d] === 0; d--)
    ;
  if (h > d && (h = d), d === 0)
    return a[r++] = 1 << 24 | 64 << 16 | 0, a[r++] = 1 << 24 | 64 << 16 | 0, o.bits = 1, 0;
  for (b = 1; b < d && g[b] === 0; b++)
    ;
  for (h < b && (h = b), C = 1, f = 1; f <= he; f++)
    if (C <<= 1, C -= g[f], C < 0)
      return -1;
  if (C > 0 && (t === Kt || d !== 1))
    return -1;
  for (N[1] = 0, f = 1; f < he; f++)
    N[f + 1] = N[f] + g[f];
  for (l = 0; l < n; l++)
    i[e + l] !== 0 && (s[N[i[e + l]]++] = l);
  if (t === Kt ? (w = te = s, m = 20) : t === ot ? (w = qa, te = Ja, m = 257) : (w = Qa, te = er, m = 0), E = 0, l = 0, f = b, _ = r, p = h, k = 0, y = -1, v = 1 << h, S = v - 1, t === ot && v > jt || t === Yt && v > Gt)
    return 1;
  for (; ; ) {
    Ot = f - k, s[l] + 1 < m ? (Be = 0, Ue = s[l]) : s[l] >= m ? (Be = te[s[l] - m], Ue = w[s[l] - m]) : (Be = 32 + 64, Ue = 0), $ = 1 << f - k, c = 1 << p, b = c;
    do
      c -= $, a[_ + (E >> k) + c] = Ot << 24 | Be << 16 | Ue | 0;
    while (c !== 0);
    for ($ = 1 << f - 1; E & $; )
      $ >>= 1;
    if ($ !== 0 ? (E &= $ - 1, E += $) : E = 0, l++, --g[f] === 0) {
      if (f === d)
        break;
      f = i[e + s[l]];
    }
    if (f > h && (E & S) !== y) {
      for (k === 0 && (k = h), _ += b, p = f - k, C = 1 << p; p + k < d && (C -= g[p + k], !(C <= 0)); )
        p++, C <<= 1;
      if (v += 1 << p, t === ot && v > jt || t === Yt && v > Gt)
        return 1;
      y = E & S, a[y] = h << 24 | p << 16 | _ - r | 0;
    }
  }
  return E !== 0 && (a[_ + E] = f - k << 24 | 64 << 16 | 0), o.bits = h, 0;
};
var $e = tr;
const ir = 0, Gi = 1, Ki = 2, {
  Z_FINISH: Xt,
  Z_BLOCK: nr,
  Z_TREES: Pe,
  Z_OK: oe,
  Z_STREAM_END: ar,
  Z_NEED_DICT: rr,
  Z_STREAM_ERROR: F,
  Z_DATA_ERROR: Yi,
  Z_MEM_ERROR: Xi,
  Z_BUF_ERROR: sr,
  Z_DEFLATED: Wt
} = Oe, et = 16180, qt = 16181, Jt = 16182, Qt = 16183, ei = 16184, ti = 16185, ii = 16186, ni = 16187, ai = 16188, ri = 16189, We = 16190, G = 16191, ft = 16192, si = 16193, lt = 16194, oi = 16195, fi = 16196, li = 16197, hi = 16198, Me = 16199, He = 16200, ui = 16201, di = 16202, ci = 16203, _i = 16204, pi = 16205, ht = 16206, bi = 16207, gi = 16208, A = 16209, Wi = 16210, qi = 16211, or = 852, fr = 592, lr = 15, hr = lr, wi = (t) => (t >>> 24 & 255) + (t >>> 8 & 65280) + ((t & 65280) << 8) + ((t & 255) << 24);
function ur() {
  this.strm = null, this.mode = 0, this.last = !1, this.wrap = 0, this.havedict = !1, this.flags = 0, this.dmax = 0, this.check = 0, this.total = 0, this.head = null, this.wbits = 0, this.wsize = 0, this.whave = 0, this.wnext = 0, this.window = null, this.hold = 0, this.bits = 0, this.length = 0, this.offset = 0, this.extra = 0, this.lencode = null, this.distcode = null, this.lenbits = 0, this.distbits = 0, this.ncode = 0, this.nlen = 0, this.ndist = 0, this.have = 0, this.next = null, this.lens = new Uint16Array(320), this.work = new Uint16Array(288), this.lendyn = null, this.distdyn = null, this.sane = 0, this.back = 0, this.was = 0;
}
const fe = (t) => {
  if (!t)
    return 1;
  const i = t.state;
  return !i || i.strm !== t || i.mode < et || i.mode > qi ? 1 : 0;
}, Ji = (t) => {
  if (fe(t))
    return F;
  const i = t.state;
  return t.total_in = t.total_out = i.total = 0, t.msg = "", i.wrap && (t.adler = i.wrap & 1), i.mode = et, i.last = 0, i.havedict = 0, i.flags = -1, i.dmax = 32768, i.head = null, i.hold = 0, i.bits = 0, i.lencode = i.lendyn = new Int32Array(or), i.distcode = i.distdyn = new Int32Array(fr), i.sane = 1, i.back = -1, oe;
}, Qi = (t) => {
  if (fe(t))
    return F;
  const i = t.state;
  return i.wsize = 0, i.whave = 0, i.wnext = 0, Ji(t);
}, en = (t, i) => {
  let e;
  if (fe(t))
    return F;
  const n = t.state;
  return i < 0 ? (e = 0, i = -i) : (e = (i >> 4) + 5, i < 48 && (i &= 15)), i && (i < 8 || i > 15) ? F : (n.window !== null && n.wbits !== i && (n.window = null), n.wrap = e, n.wbits = i, Qi(t));
}, tn = (t, i) => {
  if (!t)
    return F;
  const e = new ur();
  t.state = e, e.strm = t, e.window = null, e.mode = et;
  const n = en(t, i);
  return n !== oe && (t.state = null), n;
}, dr = (t) => tn(t, hr);
let mi = !0, ut, dt;
const cr = (t) => {
  if (mi) {
    ut = new Int32Array(512), dt = new Int32Array(32);
    let i = 0;
    for (; i < 144; )
      t.lens[i++] = 8;
    for (; i < 256; )
      t.lens[i++] = 9;
    for (; i < 280; )
      t.lens[i++] = 7;
    for (; i < 288; )
      t.lens[i++] = 8;
    for ($e(Gi, t.lens, 0, 288, ut, 0, t.work, { bits: 9 }), i = 0; i < 32; )
      t.lens[i++] = 5;
    $e(Ki, t.lens, 0, 32, dt, 0, t.work, { bits: 5 }), mi = !1;
  }
  t.lencode = ut, t.lenbits = 9, t.distcode = dt, t.distbits = 5;
}, nn = (t, i, e, n) => {
  let a;
  const r = t.state;
  return r.window === null && (r.wsize = 1 << r.wbits, r.wnext = 0, r.whave = 0, r.window = new Uint8Array(r.wsize)), n >= r.wsize ? (r.window.set(i.subarray(e - r.wsize, e), 0), r.wnext = 0, r.whave = r.wsize) : (a = r.wsize - r.wnext, a > n && (a = n), r.window.set(i.subarray(e - n, e - n + a), r.wnext), n -= a, n ? (r.window.set(i.subarray(e - n, e), 0), r.wnext = n, r.whave = r.wsize) : (r.wnext += a, r.wnext === r.wsize && (r.wnext = 0), r.whave < r.wsize && (r.whave += a))), 0;
}, _r = (t, i) => {
  let e, n, a, r, s, o, u, f, l, b, d, h, p, k, C = 0, v, E, $, c, y, S, _, w;
  const m = new Uint8Array(4);
  let g, N;
  const te = (
    /* permutation of code lengths */
    new Uint8Array([16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15])
  );
  if (fe(t) || !t.output || !t.input && t.avail_in !== 0)
    return F;
  e = t.state, e.mode === G && (e.mode = ft), s = t.next_out, a = t.output, u = t.avail_out, r = t.next_in, n = t.input, o = t.avail_in, f = e.hold, l = e.bits, b = o, d = u, w = oe;
  e:
    for (; ; )
      switch (e.mode) {
        case et:
          if (e.wrap === 0) {
            e.mode = ft;
            break;
          }
          for (; l < 16; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          if (e.wrap & 2 && f === 35615) {
            e.wbits === 0 && (e.wbits = 15), e.check = 0, m[0] = f & 255, m[1] = f >>> 8 & 255, e.check = z(e.check, m, 2, 0), f = 0, l = 0, e.mode = qt;
            break;
          }
          if (e.head && (e.head.done = !1), !(e.wrap & 1) || /* check if zlib header allowed */
          (((f & 255) << 8) + (f >> 8)) % 31) {
            t.msg = "incorrect header check", e.mode = A;
            break;
          }
          if ((f & 15) !== Wt) {
            t.msg = "unknown compression method", e.mode = A;
            break;
          }
          if (f >>>= 4, l -= 4, _ = (f & 15) + 8, e.wbits === 0 && (e.wbits = _), _ > 15 || _ > e.wbits) {
            t.msg = "invalid window size", e.mode = A;
            break;
          }
          e.dmax = 1 << e.wbits, e.flags = 0, t.adler = e.check = 1, e.mode = f & 512 ? ri : G, f = 0, l = 0;
          break;
        case qt:
          for (; l < 16; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          if (e.flags = f, (e.flags & 255) !== Wt) {
            t.msg = "unknown compression method", e.mode = A;
            break;
          }
          if (e.flags & 57344) {
            t.msg = "unknown header flags set", e.mode = A;
            break;
          }
          e.head && (e.head.text = f >> 8 & 1), e.flags & 512 && e.wrap & 4 && (m[0] = f & 255, m[1] = f >>> 8 & 255, e.check = z(e.check, m, 2, 0)), f = 0, l = 0, e.mode = Jt;
        case Jt:
          for (; l < 32; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          e.head && (e.head.time = f), e.flags & 512 && e.wrap & 4 && (m[0] = f & 255, m[1] = f >>> 8 & 255, m[2] = f >>> 16 & 255, m[3] = f >>> 24 & 255, e.check = z(e.check, m, 4, 0)), f = 0, l = 0, e.mode = Qt;
        case Qt:
          for (; l < 16; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          e.head && (e.head.xflags = f & 255, e.head.os = f >> 8), e.flags & 512 && e.wrap & 4 && (m[0] = f & 255, m[1] = f >>> 8 & 255, e.check = z(e.check, m, 2, 0)), f = 0, l = 0, e.mode = ei;
        case ei:
          if (e.flags & 1024) {
            for (; l < 16; ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            e.length = f, e.head && (e.head.extra_len = f), e.flags & 512 && e.wrap & 4 && (m[0] = f & 255, m[1] = f >>> 8 & 255, e.check = z(e.check, m, 2, 0)), f = 0, l = 0;
          } else
            e.head && (e.head.extra = null);
          e.mode = ti;
        case ti:
          if (e.flags & 1024 && (h = e.length, h > o && (h = o), h && (e.head && (_ = e.head.extra_len - e.length, e.head.extra || (e.head.extra = new Uint8Array(e.head.extra_len)), e.head.extra.set(
            n.subarray(
              r,
              // extra field is limited to 65536 bytes
              // - no need for additional size check
              r + h
            ),
            /*len + copy > state.head.extra_max - len ? state.head.extra_max : copy,*/
            _
          )), e.flags & 512 && e.wrap & 4 && (e.check = z(e.check, n, h, r)), o -= h, r += h, e.length -= h), e.length))
            break e;
          e.length = 0, e.mode = ii;
        case ii:
          if (e.flags & 2048) {
            if (o === 0)
              break e;
            h = 0;
            do
              _ = n[r + h++], e.head && _ && e.length < 65536 && (e.head.name += String.fromCharCode(_));
            while (_ && h < o);
            if (e.flags & 512 && e.wrap & 4 && (e.check = z(e.check, n, h, r)), o -= h, r += h, _)
              break e;
          } else
            e.head && (e.head.name = null);
          e.length = 0, e.mode = ni;
        case ni:
          if (e.flags & 4096) {
            if (o === 0)
              break e;
            h = 0;
            do
              _ = n[r + h++], e.head && _ && e.length < 65536 && (e.head.comment += String.fromCharCode(_));
            while (_ && h < o);
            if (e.flags & 512 && e.wrap & 4 && (e.check = z(e.check, n, h, r)), o -= h, r += h, _)
              break e;
          } else
            e.head && (e.head.comment = null);
          e.mode = ai;
        case ai:
          if (e.flags & 512) {
            for (; l < 16; ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            if (e.wrap & 4 && f !== (e.check & 65535)) {
              t.msg = "header crc mismatch", e.mode = A;
              break;
            }
            f = 0, l = 0;
          }
          e.head && (e.head.hcrc = e.flags >> 9 & 1, e.head.done = !0), t.adler = e.check = 0, e.mode = G;
          break;
        case ri:
          for (; l < 32; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          t.adler = e.check = wi(f), f = 0, l = 0, e.mode = We;
        case We:
          if (e.havedict === 0)
            return t.next_out = s, t.avail_out = u, t.next_in = r, t.avail_in = o, e.hold = f, e.bits = l, rr;
          t.adler = e.check = 1, e.mode = G;
        case G:
          if (i === nr || i === Pe)
            break e;
        case ft:
          if (e.last) {
            f >>>= l & 7, l -= l & 7, e.mode = ht;
            break;
          }
          for (; l < 3; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          switch (e.last = f & 1, f >>>= 1, l -= 1, f & 3) {
            case 0:
              e.mode = si;
              break;
            case 1:
              if (cr(e), e.mode = Me, i === Pe) {
                f >>>= 2, l -= 2;
                break e;
              }
              break;
            case 2:
              e.mode = fi;
              break;
            case 3:
              t.msg = "invalid block type", e.mode = A;
          }
          f >>>= 2, l -= 2;
          break;
        case si:
          for (f >>>= l & 7, l -= l & 7; l < 32; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          if ((f & 65535) !== (f >>> 16 ^ 65535)) {
            t.msg = "invalid stored block lengths", e.mode = A;
            break;
          }
          if (e.length = f & 65535, f = 0, l = 0, e.mode = lt, i === Pe)
            break e;
        case lt:
          e.mode = oi;
        case oi:
          if (h = e.length, h) {
            if (h > o && (h = o), h > u && (h = u), h === 0)
              break e;
            a.set(n.subarray(r, r + h), s), o -= h, r += h, u -= h, s += h, e.length -= h;
            break;
          }
          e.mode = G;
          break;
        case fi:
          for (; l < 14; ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          if (e.nlen = (f & 31) + 257, f >>>= 5, l -= 5, e.ndist = (f & 31) + 1, f >>>= 5, l -= 5, e.ncode = (f & 15) + 4, f >>>= 4, l -= 4, e.nlen > 286 || e.ndist > 30) {
            t.msg = "too many length or distance symbols", e.mode = A;
            break;
          }
          e.have = 0, e.mode = li;
        case li:
          for (; e.have < e.ncode; ) {
            for (; l < 3; ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            e.lens[te[e.have++]] = f & 7, f >>>= 3, l -= 3;
          }
          for (; e.have < 19; )
            e.lens[te[e.have++]] = 0;
          if (e.lencode = e.lendyn, e.lenbits = 7, g = { bits: e.lenbits }, w = $e(ir, e.lens, 0, 19, e.lencode, 0, e.work, g), e.lenbits = g.bits, w) {
            t.msg = "invalid code lengths set", e.mode = A;
            break;
          }
          e.have = 0, e.mode = hi;
        case hi:
          for (; e.have < e.nlen + e.ndist; ) {
            for (; C = e.lencode[f & (1 << e.lenbits) - 1], v = C >>> 24, E = C >>> 16 & 255, $ = C & 65535, !(v <= l); ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            if ($ < 16)
              f >>>= v, l -= v, e.lens[e.have++] = $;
            else {
              if ($ === 16) {
                for (N = v + 2; l < N; ) {
                  if (o === 0)
                    break e;
                  o--, f += n[r++] << l, l += 8;
                }
                if (f >>>= v, l -= v, e.have === 0) {
                  t.msg = "invalid bit length repeat", e.mode = A;
                  break;
                }
                _ = e.lens[e.have - 1], h = 3 + (f & 3), f >>>= 2, l -= 2;
              } else if ($ === 17) {
                for (N = v + 3; l < N; ) {
                  if (o === 0)
                    break e;
                  o--, f += n[r++] << l, l += 8;
                }
                f >>>= v, l -= v, _ = 0, h = 3 + (f & 7), f >>>= 3, l -= 3;
              } else {
                for (N = v + 7; l < N; ) {
                  if (o === 0)
                    break e;
                  o--, f += n[r++] << l, l += 8;
                }
                f >>>= v, l -= v, _ = 0, h = 11 + (f & 127), f >>>= 7, l -= 7;
              }
              if (e.have + h > e.nlen + e.ndist) {
                t.msg = "invalid bit length repeat", e.mode = A;
                break;
              }
              for (; h--; )
                e.lens[e.have++] = _;
            }
          }
          if (e.mode === A)
            break;
          if (e.lens[256] === 0) {
            t.msg = "invalid code -- missing end-of-block", e.mode = A;
            break;
          }
          if (e.lenbits = 9, g = { bits: e.lenbits }, w = $e(Gi, e.lens, 0, e.nlen, e.lencode, 0, e.work, g), e.lenbits = g.bits, w) {
            t.msg = "invalid literal/lengths set", e.mode = A;
            break;
          }
          if (e.distbits = 6, e.distcode = e.distdyn, g = { bits: e.distbits }, w = $e(Ki, e.lens, e.nlen, e.ndist, e.distcode, 0, e.work, g), e.distbits = g.bits, w) {
            t.msg = "invalid distances set", e.mode = A;
            break;
          }
          if (e.mode = Me, i === Pe)
            break e;
        case Me:
          e.mode = He;
        case He:
          if (o >= 6 && u >= 258) {
            t.next_out = s, t.avail_out = u, t.next_in = r, t.avail_in = o, e.hold = f, e.bits = l, Wa(t, d), s = t.next_out, a = t.output, u = t.avail_out, r = t.next_in, n = t.input, o = t.avail_in, f = e.hold, l = e.bits, e.mode === G && (e.back = -1);
            break;
          }
          for (e.back = 0; C = e.lencode[f & (1 << e.lenbits) - 1], v = C >>> 24, E = C >>> 16 & 255, $ = C & 65535, !(v <= l); ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          if (E && !(E & 240)) {
            for (c = v, y = E, S = $; C = e.lencode[S + ((f & (1 << c + y) - 1) >> c)], v = C >>> 24, E = C >>> 16 & 255, $ = C & 65535, !(c + v <= l); ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            f >>>= c, l -= c, e.back += c;
          }
          if (f >>>= v, l -= v, e.back += v, e.length = $, E === 0) {
            e.mode = pi;
            break;
          }
          if (E & 32) {
            e.back = -1, e.mode = G;
            break;
          }
          if (E & 64) {
            t.msg = "invalid literal/length code", e.mode = A;
            break;
          }
          e.extra = E & 15, e.mode = ui;
        case ui:
          if (e.extra) {
            for (N = e.extra; l < N; ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            e.length += f & (1 << e.extra) - 1, f >>>= e.extra, l -= e.extra, e.back += e.extra;
          }
          e.was = e.length, e.mode = di;
        case di:
          for (; C = e.distcode[f & (1 << e.distbits) - 1], v = C >>> 24, E = C >>> 16 & 255, $ = C & 65535, !(v <= l); ) {
            if (o === 0)
              break e;
            o--, f += n[r++] << l, l += 8;
          }
          if (!(E & 240)) {
            for (c = v, y = E, S = $; C = e.distcode[S + ((f & (1 << c + y) - 1) >> c)], v = C >>> 24, E = C >>> 16 & 255, $ = C & 65535, !(c + v <= l); ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            f >>>= c, l -= c, e.back += c;
          }
          if (f >>>= v, l -= v, e.back += v, E & 64) {
            t.msg = "invalid distance code", e.mode = A;
            break;
          }
          e.offset = $, e.extra = E & 15, e.mode = ci;
        case ci:
          if (e.extra) {
            for (N = e.extra; l < N; ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            e.offset += f & (1 << e.extra) - 1, f >>>= e.extra, l -= e.extra, e.back += e.extra;
          }
          if (e.offset > e.dmax) {
            t.msg = "invalid distance too far back", e.mode = A;
            break;
          }
          e.mode = _i;
        case _i:
          if (u === 0)
            break e;
          if (h = d - u, e.offset > h) {
            if (h = e.offset - h, h > e.whave && e.sane) {
              t.msg = "invalid distance too far back", e.mode = A;
              break;
            }
            h > e.wnext ? (h -= e.wnext, p = e.wsize - h) : p = e.wnext - h, h > e.length && (h = e.length), k = e.window;
          } else
            k = a, p = s - e.offset, h = e.length;
          h > u && (h = u), u -= h, e.length -= h;
          do
            a[s++] = k[p++];
          while (--h);
          e.length === 0 && (e.mode = He);
          break;
        case pi:
          if (u === 0)
            break e;
          a[s++] = e.length, u--, e.mode = He;
          break;
        case ht:
          if (e.wrap) {
            for (; l < 32; ) {
              if (o === 0)
                break e;
              o--, f |= n[r++] << l, l += 8;
            }
            if (d -= u, t.total_out += d, e.total += d, e.wrap & 4 && d && (t.adler = e.check = /*UPDATE_CHECK(state.check, put - _out, _out);*/
            e.flags ? z(e.check, a, d, s - d) : Te(e.check, a, d, s - d)), d = u, e.wrap & 4 && (e.flags ? f : wi(f)) !== e.check) {
              t.msg = "incorrect data check", e.mode = A;
              break;
            }
            f = 0, l = 0;
          }
          e.mode = bi;
        case bi:
          if (e.wrap && e.flags) {
            for (; l < 32; ) {
              if (o === 0)
                break e;
              o--, f += n[r++] << l, l += 8;
            }
            if (e.wrap & 4 && f !== (e.total & 4294967295)) {
              t.msg = "incorrect length check", e.mode = A;
              break;
            }
            f = 0, l = 0;
          }
          e.mode = gi;
        case gi:
          w = ar;
          break e;
        case A:
          w = Yi;
          break e;
        case Wi:
          return Xi;
        case qi:
        default:
          return F;
      }
  return t.next_out = s, t.avail_out = u, t.next_in = r, t.avail_in = o, e.hold = f, e.bits = l, (e.wsize || d !== t.avail_out && e.mode < A && (e.mode < ht || i !== Xt)) && nn(t, t.output, t.next_out, d - t.avail_out), b -= t.avail_in, d -= t.avail_out, t.total_in += b, t.total_out += d, e.total += d, e.wrap & 4 && d && (t.adler = e.check = /*UPDATE_CHECK(state.check, strm.next_out - _out, _out);*/
  e.flags ? z(e.check, a, d, t.next_out - d) : Te(e.check, a, d, t.next_out - d)), t.data_type = e.bits + (e.last ? 64 : 0) + (e.mode === G ? 128 : 0) + (e.mode === Me || e.mode === lt ? 256 : 0), (b === 0 && d === 0 || i === Xt) && w === oe && (w = sr), w;
}, pr = (t) => {
  if (fe(t))
    return F;
  let i = t.state;
  return i.window && (i.window = null), t.state = null, oe;
}, br = (t, i) => {
  if (fe(t))
    return F;
  const e = t.state;
  return e.wrap & 2 ? (e.head = i, i.done = !1, oe) : F;
}, gr = (t, i) => {
  const e = i.length;
  let n, a, r;
  return fe(t) || (n = t.state, n.wrap !== 0 && n.mode !== We) ? F : n.mode === We && (a = 1, a = Te(a, i, e, 0), a !== n.check) ? Yi : (r = nn(t, i, e, e), r ? (n.mode = Wi, Xi) : (n.havedict = 1, oe));
};
var wr = Qi, mr = en, yr = Ji, vr = dr, kr = tn, $r = _r, Cr = pr, Er = br, Nr = gr, Sr = "pako inflate (from Nodeca project)", Y = {
  inflateReset: wr,
  inflateReset2: mr,
  inflateResetKeep: yr,
  inflateInit: vr,
  inflateInit2: kr,
  inflate: $r,
  inflateEnd: Cr,
  inflateGetHeader: Er,
  inflateSetDictionary: Nr,
  inflateInfo: Sr
};
function Tr() {
  this.text = 0, this.time = 0, this.xflags = 0, this.os = 0, this.extra = null, this.extra_len = 0, this.name = "", this.comment = "", this.hcrc = 0, this.done = !1;
}
var xr = Tr;
const an = Object.prototype.toString, {
  Z_NO_FLUSH: Ir,
  Z_FINISH: Ar,
  Z_OK: Ae,
  Z_STREAM_END: ct,
  Z_NEED_DICT: _t,
  Z_STREAM_ERROR: zr,
  Z_DATA_ERROR: yi,
  Z_MEM_ERROR: Or
} = Oe;
function De(t) {
  this.options = Qe.assign({
    chunkSize: 1024 * 64,
    windowBits: 15,
    to: ""
  }, t || {});
  const i = this.options;
  i.raw && i.windowBits >= 0 && i.windowBits < 16 && (i.windowBits = -i.windowBits, i.windowBits === 0 && (i.windowBits = -15)), i.windowBits >= 0 && i.windowBits < 16 && !(t && t.windowBits) && (i.windowBits += 32), i.windowBits > 15 && i.windowBits < 48 && (i.windowBits & 15 || (i.windowBits |= 15)), this.err = 0, this.msg = "", this.ended = !1, this.chunks = [], this.strm = new Hi(), this.strm.avail_out = 0;
  let e = Y.inflateInit2(
    this.strm,
    i.windowBits
  );
  if (e !== Ae)
    throw new Error(de[e]);
  if (this.header = new xr(), Y.inflateGetHeader(this.strm, this.header), i.dictionary && (typeof i.dictionary == "string" ? i.dictionary = Ie.string2buf(i.dictionary) : an.call(i.dictionary) === "[object ArrayBuffer]" && (i.dictionary = new Uint8Array(i.dictionary)), i.raw && (e = Y.inflateSetDictionary(this.strm, i.dictionary), e !== Ae)))
    throw new Error(de[e]);
}
De.prototype.push = function(t, i) {
  const e = this.strm, n = this.options.chunkSize, a = this.options.dictionary;
  let r, s, o;
  if (this.ended)
    return !1;
  for (i === ~~i ? s = i : s = i === !0 ? Ar : Ir, an.call(t) === "[object ArrayBuffer]" ? e.input = new Uint8Array(t) : e.input = t, e.next_in = 0, e.avail_in = e.input.length; ; ) {
    for (e.avail_out === 0 && (e.output = new Uint8Array(n), e.next_out = 0, e.avail_out = n), r = Y.inflate(e, s), r === _t && a && (r = Y.inflateSetDictionary(e, a), r === Ae ? r = Y.inflate(e, s) : r === yi && (r = _t)); e.avail_in > 0 && r === ct && e.state.wrap > 0 && t[e.next_in] !== 0; )
      Y.inflateReset(e), r = Y.inflate(e, s);
    switch (r) {
      case zr:
      case yi:
      case _t:
      case Or:
        return this.onEnd(r), this.ended = !0, !1;
    }
    if (o = e.avail_out, e.next_out && (e.avail_out === 0 || r === ct))
      if (this.options.to === "string") {
        let u = Ie.utf8border(e.output, e.next_out), f = e.next_out - u, l = Ie.buf2string(e.output, u);
        e.next_out = f, e.avail_out = n - f, f && e.output.set(e.output.subarray(u, u + f), 0), this.onData(l);
      } else
        this.onData(e.output.length === e.next_out ? e.output : e.output.subarray(0, e.next_out));
    if (!(r === Ae && o === 0)) {
      if (r === ct)
        return r = Y.inflateEnd(this.strm), this.onEnd(r), this.ended = !0, !0;
      if (e.avail_in === 0)
        break;
    }
  }
  return !0;
};
De.prototype.onData = function(t) {
  this.chunks.push(t);
};
De.prototype.onEnd = function(t) {
  t === Ae && (this.options.to === "string" ? this.result = this.chunks.join("") : this.result = Qe.flattenChunks(this.chunks)), this.chunks = [], this.err = t, this.msg = this.strm.msg;
};
function zt(t, i) {
  const e = new De(i);
  if (e.push(t), e.err)
    throw e.msg || de[e.err];
  return e.result;
}
function Rr(t, i) {
  return i = i || {}, i.raw = !0, zt(t, i);
}
var Dr = De, Br = zt, Ur = Rr, Vr = zt, Fr = Oe, Zr = {
  Inflate: Dr,
  inflate: Br,
  inflateRaw: Ur,
  ungzip: Vr,
  constants: Fr
};
const { Inflate: es, inflate: ts, inflateRaw: Lr, ungzip: is } = Zr;
var Pr = Lr;
function Mr(t) {
  return Pr(t.subarray(2));
}
class Hr extends Error {
  constructor(i) {
    super(i), this.code = "ERR_ABORTED";
  }
}
function jr(t) {
  t.sort((a, r) => Number(a.offset) - Number(r.offset));
  const i = [];
  let e, n;
  for (let a = 0; a < t.length; a += 1)
    e && n && Number(t[a].offset) - n <= 2e3 ? (e.length = BigInt(Number(e.length) + Number(t[a].length) - n + Number(t[a].offset)), e.blocks.push(t[a])) : i.push(e = {
      blocks: [t[a]],
      length: t[a].length,
      offset: t[a].offset
    }), n = Number(e.offset) + Number(e.length);
  return i;
}
function je(t) {
  if (t && t.aborted) {
    if (typeof DOMException < "u")
      throw new DOMException("aborted", "AbortError");
    {
      const i = new Hr("aborted");
      throw i.code = "ERR_ABORTED", i;
    }
  }
}
const Gr = 1, Kr = 2, Yr = 3;
function pt(t, i, e, n) {
  return t < n && i >= e;
}
function Xr(t) {
  const i = t ? "big" : "little", e = new x().endianess(i).uint32("chromId").uint32("start").uint32("end").uint32("validCnt").floatle("minScore").floatle("maxScore").floatle("sumData").floatle("sumSqData").saveOffset("offset"), n = new x().endianess(i).uint8("isLeaf").skip(1).uint16("cnt").choice({
    tag: "isLeaf",
    choices: {
      1: new x().endianess(i).array("blocksToFetch", {
        length: "cnt",
        type: new x().endianess(i).uint32("startChrom").uint32("startBase").uint32("endChrom").uint32("endBase").uint64("blockOffset").uint64("blockSize").saveOffset("offset")
      }),
      0: new x().array("recurOffsets", {
        length: "cnt",
        type: new x().endianess(i).uint32("startChrom").uint32("startBase").uint32("endChrom").uint32("endBase").uint64("blockOffset").saveOffset("offset")
      })
    }
  }), a = new x().endianess(i).uint32("chromId").int32("start").int32("end").string("rest", {
    zeroTerminated: !0
  }).saveOffset("offset");
  return {
    bigWigParser: new x().endianess(i).skip(4).int32("blockStart").skip(4).uint32("itemStep").uint32("itemSpan").uint8("blockType").skip(1).uint16("itemCount").choice({
      tag: "blockType",
      choices: {
        [Yr]: new x().array("items", {
          length: "itemCount",
          type: new x().floatle("score")
        }),
        [Kr]: new x().array("items", {
          length: "itemCount",
          type: new x().endianess(i).int32("start").floatle("score")
        }),
        [Gr]: new x().array("items", {
          length: "itemCount",
          type: new x().endianess(i).int32("start").int32("end").floatle("score")
        })
      }
    }),
    bigBedParser: a,
    summaryParser: e,
    leafParser: n
  };
}
class Wr {
  constructor(i, e, n, a, r, s) {
    if (this.bbi = i, this.refsByName = e, this.cirTreeOffset = n, this.isBigEndian = a, this.isCompressed = r, this.blockType = s, this.featureCache = new rn({
      cache: new ln({ maxSize: 1e3 }),
      fill: async (u, f) => {
        const l = Number(u.length), b = Number(u.offset), { buffer: d } = await this.bbi.read(Ke.Buffer.alloc(l), 0, l, b, {
          signal: f
        });
        return d;
      }
    }), !(n >= 0))
      throw new Error("invalid cirTreeOffset!");
    const o = Xr(a);
    this.leafParser = o.leafParser, this.bigBedParser = o.bigBedParser;
  }
  async readWigData(i, e, n, a, r) {
    try {
      const { refsByName: s, bbi: o, cirTreeOffset: u, isBigEndian: f } = this, l = s[i];
      l === void 0 && a.complete();
      const b = { chrId: l, start: e, end: n };
      this.cirTreePromise || (this.cirTreePromise = o.read(Ke.Buffer.alloc(48), 0, 48, Number(u), r));
      const { buffer: d } = await this.cirTreePromise, h = f ? d.readUInt32BE(4) : d.readUInt32LE(4);
      let p = [], k = 0;
      const C = (c, y, S) => {
        try {
          const _ = c.subarray(y), w = this.leafParser.parse(_);
          if (w.blocksToFetch && (p = p.concat(w.blocksToFetch.filter(v).map((m) => ({
            offset: m.blockOffset,
            length: m.blockSize
          })))), w.recurOffsets) {
            const m = w.recurOffsets.filter(v).map((g) => Number(g.blockOffset));
            m.length > 0 && $(m, S + 1);
          }
        } catch (_) {
          a.error(_);
        }
      }, v = (c) => {
        const { startChrom: y, startBase: S, endChrom: _, endBase: w } = c;
        return (y < l || y === l && S <= n) && (_ > l || _ === l && w >= e);
      }, E = async (c, y, S) => {
        try {
          const _ = y.max() - y.min(), w = y.min(), m = await this.featureCache.get(`${_}_${w}`, { length: _, offset: w }, r.signal);
          for (let g = 0; g < c.length; g += 1)
            y.contains(c[g]) && (C(m, c[g] - w, S), k -= 1, k === 0 && this.readFeatures(a, p, { ...r, request: b }));
        } catch (_) {
          a.error(_);
        }
      }, $ = (c, y) => {
        try {
          k += c.length;
          const S = 4 + Number(h) * 32;
          let _ = new X(c[0], c[0] + S);
          for (let w = 1; w < c.length; w += 1) {
            const m = new X(c[w], c[w] + S);
            _ = _.union(m);
          }
          _.getRanges().map((w) => E(c, w, y));
        } catch (S) {
          a.error(S);
        }
      };
      return $([Number(u) + 48], 1);
    } catch (s) {
      a.error(s);
    }
  }
  parseSummaryBlock(i, e, n) {
    const a = [];
    let r = e;
    const s = new DataView(i.buffer, i.byteOffset, i.length);
    for (; r < i.byteLength; ) {
      const o = s.getUint32(r, !0);
      r += 4;
      const u = s.getUint32(r, !0);
      r += 4;
      const f = s.getUint32(r, !0);
      r += 4;
      const l = s.getUint32(r, !0);
      r += 4;
      const b = s.getFloat32(r, !0);
      r += 4;
      const d = s.getFloat32(r, !0);
      r += 4;
      const h = s.getFloat32(r, !0);
      r += 4, r += 4, (!n || o === n.chrId && pt(u, f, n.start, n.end)) && a.push({
        start: u,
        end: f,
        maxScore: d,
        minScore: b,
        summary: !0,
        score: h / (l || 1)
      });
    }
    return a;
  }
  parseBigBedBlock(i, e, n, a) {
    const r = [];
    let s = e;
    for (; s < i.byteLength; ) {
      const o = this.bigBedParser.parse(i.subarray(s));
      r.push({ ...o, uniqueId: `bb-${n + s}` }), s += o.offset;
    }
    return a ? r.filter((o) => pt(o.start, o.end, a.start, a.end)) : r;
  }
  parseBigWigBlock(i, e, n) {
    const a = i.subarray(e), r = new DataView(a.buffer, a.byteOffset, a.length);
    let s = 0;
    s += 4;
    const o = r.getInt32(s, !0);
    s += 8;
    const u = r.getUint32(s, !0);
    s += 4;
    const f = r.getUint32(s, !0);
    s += 4;
    const l = r.getUint8(s);
    s += 2;
    const b = r.getUint16(s, !0);
    s += 2;
    const d = new Array(b);
    switch (l) {
      case 1:
        for (let h = 0; h < b; h++) {
          const p = r.getInt32(s, !0);
          s += 4;
          const k = r.getInt32(s, !0);
          s += 4;
          const C = r.getFloat32(s, !0);
          s += 4, d[h] = { start: p, end: k, score: C };
        }
        break;
      case 2:
        for (let h = 0; h < b; h++) {
          const p = r.getInt32(s, !0);
          s += 4;
          const k = r.getFloat32(s, !0);
          s += 4, d[h] = { score: k, start: p, end: p + f };
        }
        break;
      case 3:
        for (let h = 0; h < b; h++) {
          const p = r.getFloat32(s, !0);
          s += 4;
          const k = o + h * u;
          d[h] = { score: p, start: k, end: k + f };
        }
        break;
    }
    return n ? d.filter((h) => pt(h.start, h.end, n.start, n.end)) : d;
  }
  async readFeatures(i, e, n = {}) {
    try {
      const { blockType: a, isCompressed: r } = this, { signal: s, request: o } = n, u = jr(e);
      je(s), await Promise.all(u.map(async (f) => {
        je(s);
        const { length: l, offset: b } = f, d = await this.featureCache.get(`${l}_${b}`, f, s);
        f.blocks.forEach((h) => {
          je(s);
          let p = Number(h.offset) - Number(f.offset), k = d;
          switch (r && (k = Mr(d.subarray(p)), p = 0), je(s), a) {
            case "summary":
              i.next(this.parseSummaryBlock(k, p, o));
              break;
            case "bigwig":
              i.next(this.parseBigWigBlock(k, p, o));
              break;
            case "bigbed":
              i.next(this.parseBigBedBlock(k, p, Number(h.offset) * 256, o));
              break;
            default:
              console.warn(`Don't know what to do with ${a}`);
          }
        });
      })), i.complete();
    } catch (a) {
      i.error(a);
    }
  }
}
const vi = -2003829722, bt = -2021002517;
function qr(t) {
  return new TextDecoder().decode(t);
}
function ki(t) {
  const i = t ? "big" : "little", e = new x().endianess(i).int32("magic").uint16("version").uint16("numZoomLevels").uint64("chromTreeOffset").uint64("unzoomedDataOffset").uint64("unzoomedIndexOffset").uint16("fieldCount").uint16("definedFieldCount").uint64("asOffset").uint64("totalSummaryOffset").uint32("uncompressBufSize").uint64("extHeaderOffset").array("zoomLevels", {
    length: "numZoomLevels",
    type: new x().endianess(i).uint32("reductionLevel").uint32("reserved").uint64("dataOffset").uint64("indexOffset")
  }), n = new x().endianess(i).uint64("basesCovered").doublele("scoreMin").doublele("scoreMax").doublele("scoreSum").doublele("scoreSumSquares"), a = new x().endianess(i).uint32("magic").uint32("blockSize").uint32("keySize").uint32("valSize").uint64("itemCount"), r = new x().endianess(i).uint8("isLeafNode").skip(1).uint16("cnt").saveOffset("offset");
  return {
    chromTreeParser: a,
    totalSummaryParser: n,
    headerParser: e,
    isLeafNode: r
  };
}
class ns {
  /* fetch and parse header information from a bigwig or bigbed file
   * @param abortSignal - abort the operation, can be null
   * @return a Header object
   */
  getHeader(i = {}) {
    const e = "aborted" in i ? { signal: i } : i;
    return this.headerP || (this.headerP = this._getHeader(e).catch((n) => {
      throw this.headerP = void 0, n;
    })), this.headerP;
  }
  /*
   * @param filehandle - a filehandle from generic-filehandle or implementing something similar to the node10 fs.promises API
   * @param path - a Local file path as a string
   * @param url - a URL string
   * @param renameRefSeqs - an optional method to rename the internal reference sequences using a mapping function
   */
  constructor(i = {}) {
    const { filehandle: e, renameRefSeqs: n = (s) => s, path: a, url: r } = i;
    if (this.renameRefSeqs = n, e)
      this.bbi = e;
    else if (r)
      this.bbi = new sn(r);
    else if (a)
      this.bbi = new on(a);
    else
      throw new Error("no file given");
  }
  async _getHeader(i) {
    const e = await this._getMainHeader(i), n = await this._readChromTree(e, i);
    return { ...e, ...n };
  }
  async _getMainHeader(i, e = 2e3) {
    const { buffer: n } = await this.bbi.read(Ke.Buffer.alloc(e), 0, e, 0, i), a = this._isBigEndian(n), r = ki(a), s = r.headerParser.parse(n), { magic: o, asOffset: u, totalSummaryOffset: f } = s;
    if (s.fileType = o === bt ? "bigbed" : "bigwig", u > e || f > e)
      return this._getMainHeader(i, e * 2);
    if (u) {
      const l = Number(s.asOffset);
      s.autoSql = qr(n.subarray(l, n.indexOf(0, l)));
    }
    if (s.totalSummaryOffset > e)
      return this._getMainHeader(i, e * 2);
    if (s.totalSummaryOffset) {
      const l = n.subarray(Number(s.totalSummaryOffset)), b = r.totalSummaryParser.parse(l);
      s.totalSummary = { ...b, basesCovered: Number(b.basesCovered) };
    }
    return { ...s, isBigEndian: a };
  }
  _isBigEndian(i) {
    let e = i.readInt32LE(0);
    if (e === vi || e === bt)
      return !1;
    if (e = i.readInt32BE(0), e === vi || e === bt)
      return !0;
    throw new Error("not a BigWig/BigBed file");
  }
  // todo: add progress if long running
  async _readChromTree(i, e) {
    const n = i.isBigEndian, a = n ? "big" : "little", r = [], s = {};
    let o = Number(i.unzoomedDataOffset);
    const u = Number(i.chromTreeOffset);
    for (; o % 4 !== 0; )
      o += 1;
    const f = o - u, { buffer: l } = await this.bbi.read(Ke.Buffer.alloc(f), 0, f, Number(u), e), b = ki(n), { keySize: d } = b.chromTreeParser.parse(l), h = new x().endianess(a).string("key", { stripNull: !0, length: d }).uint32("refId").uint32("refSize").saveOffset("offset"), p = new x().endianess(a).skip(d).uint64("childOffset").saveOffset("offset"), k = 32, C = async (v) => {
      let E = v;
      if (E >= l.length)
        throw new Error("reading beyond end of buffer");
      const $ = b.isLeafNode.parse(l.subarray(E)), { isLeafNode: c, cnt: y } = $;
      if (E += $.offset, c)
        for (let S = 0; S < y; S += 1) {
          const _ = h.parse(l.subarray(E));
          E += _.offset;
          const { key: w, refId: m, refSize: g } = _, N = { name: w, id: m, length: g };
          s[this.renameRefSeqs(w)] = m, r[m] = N;
        }
      else {
        const S = [];
        for (let _ = 0; _ < y; _ += 1) {
          const w = p.parse(l.subarray(E)), { childOffset: m } = w;
          E += w.offset, S.push(C(Number(m) - Number(u)));
        }
        await Promise.all(S);
      }
    };
    return await C(k), {
      refsByName: s,
      refsByNumber: r
    };
  }
  /*
   * fetches the "unzoomed" view of the bigwig data. this is the default for bigbed
   * @param abortSignal - a signal to optionally abort this operation
   */
  async getUnzoomedView(i) {
    const { unzoomedIndexOffset: e, refsByName: n, uncompressBufSize: a, isBigEndian: r, fileType: s } = await this.getHeader(i);
    return new Wr(this.bbi, n, e, r, a > 0, s);
  }
  /**
   * Gets features from a BigWig file
   *
   * @param refName - The chromosome name
   * @param start - The start of a region
   * @param end - The end of a region
   * @param opts - An object containing basesPerSpan (e.g. pixels per basepair) or scale used to infer the zoomLevel to use
   */
  async getFeatureStream(i, e, n, a = {
    scale: 1
  }) {
    await this.getHeader(a);
    const r = this.renameRefSeqs(i);
    let s;
    if (a.basesPerSpan ? s = await this.getView(1 / a.basesPerSpan, a) : a.scale ? s = await this.getView(a.scale, a) : s = await this.getView(1, a), !s)
      throw new Error("unable to get block view for data");
    return new Rt.Observable((o) => {
      s.readWigData(r, e, n, o, a);
    });
  }
  async getFeatures(i, e, n, a = {
    scale: 1
  }) {
    const r = await this.getFeatureStream(i, e, n, a);
    return (await Rt.firstValueFrom(r.pipe(fn.toArray()))).flat();
  }
}
export {
  ns as B,
  x as P,
  Wr as a
};
