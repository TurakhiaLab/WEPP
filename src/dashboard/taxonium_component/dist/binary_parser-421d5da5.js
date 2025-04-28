import { b3 as buffer } from "./index-8ee3baa4.js";
import { l as long } from "./long-a46ae85c.js";
var vmBrowserify = {};
(function(exports) {
  var indexOf = function(e, t) {
    if (e.indexOf)
      return e.indexOf(t);
    for (var r = 0; r < e.length; r++)
      if (e[r] === t)
        return r;
    return -1;
  }, Object_keys = function(e) {
    if (Object.keys)
      return Object.keys(e);
    var t = [];
    for (var r in e)
      t.push(r);
    return t;
  }, forEach = function(e, t) {
    if (e.forEach)
      return e.forEach(t);
    for (var r = 0; r < e.length; r++)
      t(e[r], r, e);
  }, defineProp = function() {
    try {
      return Object.defineProperty({}, "_", {}), function(e, t, r) {
        Object.defineProperty(e, t, {
          writable: !0,
          enumerable: !1,
          configurable: !0,
          value: r
        });
      };
    } catch {
      return function(t, r, o) {
        t[r] = o;
      };
    }
  }(), globals = [
    "Array",
    "Boolean",
    "Date",
    "Error",
    "EvalError",
    "Function",
    "Infinity",
    "JSON",
    "Math",
    "NaN",
    "Number",
    "Object",
    "RangeError",
    "ReferenceError",
    "RegExp",
    "String",
    "SyntaxError",
    "TypeError",
    "URIError",
    "decodeURI",
    "decodeURIComponent",
    "encodeURI",
    "encodeURIComponent",
    "escape",
    "eval",
    "isFinite",
    "isNaN",
    "parseFloat",
    "parseInt",
    "undefined",
    "unescape"
  ];
  function Context() {
  }
  Context.prototype = {};
  var Script = exports.Script = function(t) {
    if (!(this instanceof Script))
      return new Script(t);
    this.code = t;
  };
  Script.prototype.runInContext = function(e) {
    if (!(e instanceof Context))
      throw new TypeError("needs a 'context' argument.");
    var t = document.createElement("iframe");
    t.style || (t.style = {}), t.style.display = "none", document.body.appendChild(t);
    var r = t.contentWindow, o = r.eval, s = r.execScript;
    !o && s && (s.call(r, "null"), o = r.eval), forEach(Object_keys(e), function(n) {
      r[n] = e[n];
    }), forEach(globals, function(n) {
      e[n] && (r[n] = e[n]);
    });
    var i = Object_keys(r), f = o.call(r, this.code);
    return forEach(Object_keys(r), function(n) {
      (n in e || indexOf(i, n) === -1) && (e[n] = r[n]);
    }), forEach(globals, function(n) {
      n in e || defineProp(e, n, r[n]);
    }), document.body.removeChild(t), f;
  }, Script.prototype.runInThisContext = function() {
    return eval(this.code);
  }, Script.prototype.runInNewContext = function(e) {
    var t = Script.createContext(e), r = this.runInContext(t);
    return e && forEach(Object_keys(t), function(o) {
      e[o] = t[o];
    }), r;
  }, forEach(Object_keys(Script.prototype), function(e) {
    exports[e] = Script[e] = function(t) {
      var r = Script(t);
      return r[e].apply(r, [].slice.call(arguments, 1));
    };
  }), exports.isContext = function(e) {
    return e instanceof Context;
  }, exports.createScript = function(e) {
    return exports.Script(e);
  }, exports.createContext = Script.createContext = function(e) {
    var t = new Context();
    return typeof e == "object" && forEach(Object_keys(e), function(r) {
      t[r] = e[r];
    }), t;
  };
})(vmBrowserify);
var context = {};
function _typeof$1(e) {
  return typeof Symbol == "function" && typeof Symbol.iterator == "symbol" ? _typeof$1 = function(r) {
    return typeof r;
  } : _typeof$1 = function(r) {
    return r && typeof Symbol == "function" && r.constructor === Symbol && r !== Symbol.prototype ? "symbol" : typeof r;
  }, _typeof$1(e);
}
var Context$1 = function e() {
  this.code = "", this.scopes = [["vars"]], this.isAsync = !1, this.bitFields = [], this.tmpVariableCount = 0, this.references = {};
};
Context$1.prototype.generateVariable = function(e) {
  var t = [];
  for (Array.prototype.push.apply(t, this.scopes[this.scopes.length - 1]); /^\$parent\./.test(e); )
    t.pop(), e = e.replace(/^\$parent\./, "");
  return e && t.push(e), t.join(".");
};
Context$1.prototype.generateOption = function(e) {
  switch (_typeof$1(e)) {
    case "number":
      return e.toString();
    case "string":
      return this.generateVariable(e);
    case "function":
      return "(".concat(e, ").call(").concat(this.generateVariable(), ", vars)");
    default:
      return;
  }
};
Context$1.prototype.generateError = function() {
  var e = Array.prototype.slice.call(arguments), t = Context$1.interpolate.apply(this, e);
  this.isAsync ? this.pushCode("return process.nextTick(function() { callback(new Error(".concat(
    t,
    "), vars); });"
  )) : this.pushCode("throw new Error(".concat(t, ");"));
};
Context$1.prototype.generateTmpVariable = function() {
  return "$tmp".concat(this.tmpVariableCount++);
};
Context$1.prototype.pushCode = function() {
  var e = Array.prototype.slice.call(arguments);
  this.code += "".concat(Context$1.interpolate.apply(this, e), `
`);
};
Context$1.prototype.pushPath = function(e) {
  e && this.scopes[this.scopes.length - 1].push(e);
};
Context$1.prototype.popPath = function(e) {
  e && this.scopes[this.scopes.length - 1].pop();
};
Context$1.prototype.pushScope = function(e) {
  this.scopes.push([e]);
};
Context$1.prototype.popScope = function() {
  this.scopes.pop();
};
Context$1.prototype.addReference = function(e) {
  this.references[e] || (this.references[e] = { resolved: !1, requested: !1 });
};
Context$1.prototype.markResolved = function(e) {
  this.references[e].resolved = !0;
};
Context$1.prototype.markRequested = function(e) {
  e.forEach(
    function(t) {
      this.references[t].requested = !0;
    }.bind(this)
  );
};
Context$1.prototype.getUnresolvedReferences = function() {
  var e = this.references;
  return Object.keys(this.references).filter(function(t) {
    return !e[t].resolved && !e[t].requested;
  });
};
Context$1.interpolate = function(e) {
  var t = /{\d+}/g, r = e.match(t), o = Array.prototype.slice.call(arguments, 1);
  return r && r.forEach(function(s) {
    var i = parseInt(s.substr(1, s.length - 2), 10);
    e = e.replace(s, o[i].toString());
  }), e;
};
context.Context = Context$1;
function _typeof(e) {
  return typeof Symbol == "function" && typeof Symbol.iterator == "symbol" ? _typeof = function(r) {
    return typeof r;
  } : _typeof = function(r) {
    return r && typeof Symbol == "function" && r.constructor === Symbol && r !== Symbol.prototype ? "symbol" : typeof r;
  }, _typeof(e);
}
var _require = buffer, Buffer = _require.Buffer, vm = vmBrowserify, Context = context.Context, Long = long;
typeof window < "u" && (window.Buffer = Buffer);
typeof self < "u" && (self.Buffer = Buffer);
var PRIMITIVE_TYPES = {
  UInt8: 1,
  UInt16LE: 2,
  UInt16BE: 2,
  UInt32LE: 4,
  UInt32BE: 4,
  Int8: 1,
  Int16LE: 2,
  Int16BE: 2,
  Int32LE: 4,
  Int32BE: 4,
  FloatLE: 4,
  FloatBE: 4,
  DoubleLE: 8,
  DoubleBE: 8,
  UInt64: 8,
  Int64: 8
}, SPECIAL_TYPES = {
  String: null,
  Buffer: null,
  Array: null,
  Skip: null,
  Choice: null,
  Nest: null,
  Bit: null,
  Itf8: null,
  Ltf8: null
}, aliasRegistry = {}, FUNCTION_PREFIX = "___parser_", BIT_RANGE = [];
(function() {
  var e;
  for (e = 1; e <= 32; e++)
    BIT_RANGE.push(e);
})();
var NAME_MAP = {};
Object.keys(PRIMITIVE_TYPES).concat(Object.keys(SPECIAL_TYPES)).forEach(function(e) {
  NAME_MAP[e.toLowerCase()] = e;
});
var Parser = function e() {
  this.varName = "", this.type = "", this.options = {}, this.next = null, this.head = null, this.compiled = null, this.endian = "le", this.constructorFn = null, this.alias = null;
};
Parser.start = function() {
  return new Parser();
};
Object.keys(PRIMITIVE_TYPES).forEach(function(e) {
  Parser.prototype[e.toLowerCase()] = function(r, o) {
    return this.setNextParser(e.toLowerCase(), r, o);
  };
  var t = e.replace(/BE|LE/, "").toLowerCase();
  t in Parser.prototype || (Parser.prototype[t] = function(r, o) {
    return this[t + this.endian](r, o);
  });
});
BIT_RANGE.forEach(function(e) {
  Parser.prototype["bit".concat(e.toString())] = function(t, r) {
    return r || (r = {}), r.length = e, this.setNextParser("bit", t, r);
  };
});
Parser.prototype.namely = function(e) {
  return aliasRegistry[e] = this, this.alias = e, this;
};
Parser.prototype.skip = function(e, t) {
  if (t && t.assert)
    throw new Error("assert option on skip is not allowed.");
  return this.setNextParser("skip", "", { length: e });
};
Parser.prototype.string = function(e, t) {
  if (!t.zeroTerminated && !t.length && !t.greedy)
    throw new Error(
      "Neither length, zeroTerminated, nor greedy is defined for string."
    );
  if ((t.zeroTerminated || t.length) && t.greedy)
    throw new Error(
      "greedy is mutually exclusive with length and zeroTerminated for string."
    );
  if (t.stripNull && !(t.length || t.greedy))
    throw new Error(
      "Length or greedy must be defined if stripNull is defined."
    );
  return t.encoding = t.encoding || "utf8", this.setNextParser("string", e, t);
};
Parser.prototype.buffer = function(e, t) {
  if (!t.length && !t.readUntil)
    throw new Error("Length nor readUntil is defined in buffer parser");
  return this.setNextParser("buffer", e, t);
};
Parser.prototype.array = function(e, t) {
  if (!t.readUntil && !t.length && !t.lengthInBytes)
    throw new Error("Length option of array is not defined.");
  if (!t.type)
    throw new Error("Type option of array is not defined.");
  if (typeof t.type == "string" && !aliasRegistry[t.type] && Object.keys(PRIMITIVE_TYPES).indexOf(NAME_MAP[t.type]) < 0)
    throw new Error('Specified primitive type "'.concat(
      t.type,
      '" is not supported.'
    ));
  return this.setNextParser("array", e, t);
};
Parser.prototype.choice = function(e, t) {
  if (arguments.length === 1 && _typeof(e) === "object" && (t = e, e = null), !t.tag)
    throw new Error("Tag option of array is not defined.");
  if (!t.choices)
    throw new Error("Choices option of array is not defined.");
  return Object.keys(t.choices).forEach(function(r) {
    if (!t.choices[r])
      throw new Error("Choice Case ".concat(r, " of ").concat(e, " is not valid."));
    if (typeof t.choices[r] == "string" && !aliasRegistry[t.choices[r]] && Object.keys(PRIMITIVE_TYPES).indexOf(NAME_MAP[t.choices[r]]) < 0)
      throw new Error('Specified primitive type "'.concat(
        t.choices[r],
        '" is not supported.'
      ));
  }, this), this.setNextParser("choice", e, t);
};
Parser.prototype.nest = function(e, t) {
  if (arguments.length === 1 && _typeof(e) === "object" && (t = e, e = null), !t.type)
    throw new Error("Type option of nest is not defined.");
  if (!(t.type instanceof Parser) && !aliasRegistry[t.type])
    throw new Error("Type option of nest must be a Parser object.");
  if (!(t.type instanceof Parser) && !e)
    throw new Error(
      "options.type must be a object if variable name is omitted."
    );
  return this.setNextParser("nest", e, t);
};
Parser.prototype.endianess = function(e) {
  switch (e.toLowerCase()) {
    case "little":
      this.endian = "le";
      break;
    case "big":
      this.endian = "be";
      break;
    default:
      throw new Error("Invalid endianess: ".concat(e));
  }
  return this;
};
Parser.prototype.create = function(e) {
  if (!(e instanceof Function))
    throw new Error("Constructor must be a Function object.");
  return this.constructorFn = e, this;
};
Parser.prototype.getCode = function() {
  var e = new Context();
  return e.pushCode("if (!Buffer.isBuffer(buffer)) {"), e.generateError('"argument buffer is not a Buffer object"'), e.pushCode("}"), this.alias ? this.addAliasedCode(e) : this.addRawCode(e), this.alias ? e.pushCode("return {0}(0)", FUNCTION_PREFIX + this.alias) : e.pushCode("return { offset: offset, result: vars };"), e.code;
};
Parser.prototype.addRawCode = function(e) {
  e.pushCode("var offset = 0;"), this.constructorFn ? e.pushCode("var vars = new constructorFn();") : e.pushCode("var vars = {};"), this.generate(e), this.resolveReferences(e), e.pushCode("return { offset: offset, result: vars };");
};
Parser.prototype.addAliasedCode = function(e) {
  return e.pushCode("function {0}(offset) {", FUNCTION_PREFIX + this.alias), this.constructorFn ? e.pushCode("var vars = new constructorFn();") : e.pushCode("var vars = {};"), this.generate(e), e.markResolved(this.alias), this.resolveReferences(e), e.pushCode("return { offset: offset, result: vars };"), e.pushCode("}"), e;
};
Parser.prototype.resolveReferences = function(e) {
  var t = e.getUnresolvedReferences();
  e.markRequested(t), t.forEach(function(r) {
    var o = aliasRegistry[r];
    o.addAliasedCode(e);
  });
};
Parser.prototype.compile = function() {
  var e = "(function(buffer, constructorFn, Long) { ".concat(this.getCode(), " })");
  this.compiled = vm.runInThisContext(e);
};
Parser.prototype.sizeOf = function() {
  var e = NaN;
  if (Object.keys(PRIMITIVE_TYPES).indexOf(this.type) >= 0)
    e = PRIMITIVE_TYPES[this.type];
  else if (this.type === "String" && typeof this.options.length == "number")
    e = this.options.length;
  else if (this.type === "Buffer" && typeof this.options.length == "number")
    e = this.options.length;
  else if (this.type === "Array" && typeof this.options.length == "number") {
    var t = NaN;
    typeof this.options.type == "string" ? t = PRIMITIVE_TYPES[NAME_MAP[this.options.type]] : this.options.type instanceof Parser && (t = this.options.type.sizeOf()), e = this.options.length * t;
  } else
    this.type === "Skip" ? e = this.options.length : this.type === "Nest" ? e = this.options.type.sizeOf() : this.type || (e = 0);
  return this.next && (e += this.next.sizeOf()), e;
};
Parser.prototype.parse = function(e) {
  return this.compiled || this.compile(), this.compiled(e, this.constructorFn, Long);
};
Parser.prototype.setNextParser = function(e, t, r) {
  var o = new Parser();
  return o.type = NAME_MAP[e], o.varName = t, o.options = r || o.options, o.endian = this.endian, this.head ? this.head.next = o : this.next = o, this.head = o, this;
};
Parser.prototype.generate = function(e) {
  this.type && (this["generate".concat(this.type)](e), this.generateAssert(e));
  var t = e.generateVariable(this.varName);
  return this.options.formatter && this.generateFormatter(e, t, this.options.formatter), this.generateNext(e);
};
Parser.prototype.generateAssert = function(e) {
  if (this.options.assert) {
    var t = e.generateVariable(this.varName);
    switch (_typeof(this.options.assert)) {
      case "function":
        e.pushCode(
          "if (!({0}).call(vars, {1})) {",
          this.options.assert,
          t
        );
        break;
      case "number":
        e.pushCode("if ({0} !== {1}) {", this.options.assert, t);
        break;
      case "string":
        e.pushCode('if ("{0}" !== {1}) {', this.options.assert, t);
        break;
      default:
        throw new Error(
          "Assert option supports only strings, numbers and assert functions."
        );
    }
    e.generateError('"Assert error: {0} is " + {0}', t), e.pushCode("}");
  }
};
Parser.prototype.generateNext = function(e) {
  return this.next && (e = this.next.generate(e)), e;
};
Object.keys(PRIMITIVE_TYPES).forEach(function(e) {
  Parser.prototype["generate".concat(e)] = function(t) {
    e === "UInt64" ? t.pushCode(
      "{0} = Long.fromBytes(buffer.slice(offset,offset+8), true, this.endian === 'le').toNumber();",
      t.generateVariable(this.varName),
      e
    ) : e === "Int64" ? t.pushCode(
      "{0} = Long.fromBytes(buffer.slice(offset,offset+8), false, this.endian === 'le').toNumber();",
      t.generateVariable(this.varName),
      e
    ) : t.pushCode(
      "{0} = buffer.read{1}(offset);",
      t.generateVariable(this.varName),
      e
    ), t.pushCode("offset += {0};", PRIMITIVE_TYPES[e]);
  };
});
Parser.prototype.generateBit = function(e) {
  var t = JSON.parse(JSON.stringify(this));
  if (t.varName = e.generateVariable(t.varName), e.bitFields.push(t), !this.next || this.next && ["Bit", "Nest"].indexOf(this.next.type) < 0) {
    var r = 0;
    e.bitFields.forEach(function(a) {
      r += a.options.length;
    });
    var o = e.generateTmpVariable();
    if (r <= 8)
      e.pushCode("var {0} = buffer.readUInt8(offset);", o), r = 8;
    else if (r <= 16)
      e.pushCode("var {0} = buffer.readUInt16BE(offset);", o), r = 16;
    else if (r <= 24) {
      var s = e.generateTmpVariable(), i = e.generateTmpVariable();
      e.pushCode("var {0} = buffer.readUInt16BE(offset);", s), e.pushCode("var {0} = buffer.readUInt8(offset + 2);", i), e.pushCode("var {2} = ({0} << 8) | {1};", s, i, o), r = 24;
    } else if (r <= 32)
      e.pushCode("var {0} = buffer.readUInt32BE(offset);", o), r = 32;
    else
      throw new Error(
        "Currently, bit field sequence longer than 4-bytes is not supported."
      );
    e.pushCode("offset += {0};", r / 8);
    var f = 0, n = this.endian === "be";
    e.bitFields.forEach(function(a) {
      e.pushCode(
        "{0} = {1} >> {2} & {3};",
        a.varName,
        o,
        n ? r - f - a.options.length : f,
        (1 << a.options.length) - 1
      ), f += a.options.length;
    }), e.bitFields = [];
  }
};
Parser.prototype.generateSkip = function(e) {
  var t = e.generateOption(this.options.length);
  e.pushCode("offset += {0};", t);
};
Parser.prototype.generateString = function(e) {
  var t = e.generateVariable(this.varName), r = e.generateTmpVariable();
  this.options.length && this.options.zeroTerminated ? (e.pushCode("var {0} = offset;", r), e.pushCode(
    "while(buffer.readUInt8(offset++) !== 0 && offset - {0}  < {1});",
    r,
    this.options.length
  ), e.pushCode(
    "{0} = buffer.toString('{1}', {2}, offset - {2} < {3} ? offset - 1 : offset);",
    t,
    this.options.encoding,
    r,
    this.options.length
  )) : this.options.length ? (e.pushCode(
    "{0} = buffer.toString('{1}', offset, offset + {2});",
    t,
    this.options.encoding,
    e.generateOption(this.options.length)
  ), e.pushCode("offset += {0};", e.generateOption(this.options.length))) : this.options.zeroTerminated ? (e.pushCode("var {0} = offset;", r), e.pushCode("while(buffer.readUInt8(offset++) !== 0);"), e.pushCode(
    "{0} = buffer.toString('{1}', {2}, offset - 1);",
    t,
    this.options.encoding,
    r
  )) : this.options.greedy && (e.pushCode("var {0} = offset;", r), e.pushCode("while(buffer.length > offset++);"), e.pushCode(
    "{0} = buffer.toString('{1}', {2}, offset);",
    t,
    this.options.encoding,
    r
  )), this.options.stripNull && e.pushCode("{0} = {0}.replace(/\\x00+$/g, '')", t);
};
Parser.prototype.generateBuffer = function(e) {
  this.options.readUntil === "eof" ? e.pushCode(
    "{0} = buffer.slice(offset);",
    e.generateVariable(this.varName)
  ) : (e.pushCode(
    "{0} = buffer.slice(offset, offset + {1});",
    e.generateVariable(this.varName),
    e.generateOption(this.options.length)
  ), e.pushCode("offset += {0};", e.generateOption(this.options.length))), this.options.clone && e.pushCode("{0} = Buffer.from({0});", e.generateVariable(this.varName));
};
Parser.prototype.generateArray = function(e) {
  var t = e.generateOption(this.options.length), r = e.generateOption(this.options.lengthInBytes), o = this.options.type, s = e.generateTmpVariable(), i = e.generateVariable(this.varName), f = e.generateTmpVariable(), n = this.options.key, a = typeof n == "string";
  if (a ? e.pushCode("{0} = {};", i) : e.pushCode("{0} = [];", i), typeof this.options.readUntil == "function" ? e.pushCode("do {") : this.options.readUntil === "eof" ? e.pushCode("for (var {0} = 0; offset < buffer.length; {0}++) {", s) : r !== void 0 ? e.pushCode(
    "for (var {0} = offset; offset - {0} < {1}; ) {",
    s,
    r
  ) : e.pushCode("for (var {0} = 0; {0} < {1}; {0}++) {", s, t), typeof o == "string")
    if (!aliasRegistry[o])
      e.pushCode("var {0} = buffer.read{1}(offset);", f, NAME_MAP[o]), e.pushCode("offset += {0};", PRIMITIVE_TYPES[NAME_MAP[o]]);
    else {
      var p = e.generateTmpVariable();
      e.pushCode("var {0} = {1}(offset);", p, FUNCTION_PREFIX + o), e.pushCode("var {0} = {1}.result; offset = {1}.offset;", f, p), o !== this.alias && e.addReference(o);
    }
  else
    o instanceof Parser && (e.pushCode("var {0} = {};", f), e.pushScope(f), o.generate(e), e.popScope());
  a ? e.pushCode("{0}[{2}.{1}] = {2};", i, n, f) : e.pushCode("{0}.push({1});", i, f), e.pushCode("}"), typeof this.options.readUntil == "function" && e.pushCode(
    " while (!({0}).call(this, {1}, buffer.slice(offset)));",
    this.options.readUntil,
    f
  );
};
Parser.prototype.generateChoiceCase = function(e, t, r) {
  if (typeof r == "string")
    if (!aliasRegistry[r])
      e.pushCode(
        "{0} = buffer.read{1}(offset);",
        e.generateVariable(this.varName),
        NAME_MAP[r]
      ), e.pushCode("offset += {0};", PRIMITIVE_TYPES[NAME_MAP[r]]);
    else {
      var o = e.generateTmpVariable();
      e.pushCode("var {0} = {1}(offset);", o, FUNCTION_PREFIX + r), e.pushCode(
        "{0} = {1}.result; offset = {1}.offset;",
        e.generateVariable(this.varName),
        o
      ), r !== this.alias && e.addReference(r);
    }
  else
    r instanceof Parser && (e.pushPath(t), r.generate(e), e.popPath(t));
};
Parser.prototype.generateChoice = function(e) {
  var t = e.generateOption(this.options.tag);
  this.varName && e.pushCode("{0} = {};", e.generateVariable(this.varName)), e.pushCode("switch({0}) {", t), Object.keys(this.options.choices).forEach(function(r) {
    var o = this.options.choices[r];
    Number.isNaN(parseInt(r, 10)) ? e.pushCode("case '{0}':", r) : e.pushCode("case {0}:", r), this.generateChoiceCase(e, this.varName, o), e.pushCode("break;");
  }, this), e.pushCode("default:"), this.options.defaultChoice ? this.generateChoiceCase(e, this.varName, this.options.defaultChoice) : e.generateError('"Met undefined tag value " + {0} + " at choice"', t), e.pushCode("}");
};
Parser.prototype.generateNest = function(e) {
  var t = e.generateVariable(this.varName);
  if (this.options.type instanceof Parser)
    this.varName && e.pushCode("{0} = {};", t), e.pushPath(this.varName), this.options.type.generate(e), e.popPath(this.varName);
  else if (aliasRegistry[this.options.type]) {
    var r = e.generateTmpVariable();
    e.pushCode(
      "var {0} = {1}(offset);",
      r,
      FUNCTION_PREFIX + this.options.type
    ), e.pushCode("{0} = {1}.result; offset = {1}.offset;", t, r), this.options.type !== this.alias && e.addReference(this.options.type);
  }
};
Parser.prototype.generateFormatter = function(e, t, r) {
  typeof r == "function" && e.pushCode("{0} = ({1}).call(this, {0});", t, r);
};
Parser.prototype.isInteger = function() {
  return !!this.type.match(/U?Int[8|16|32][BE|LE]?|Bit\d+/);
};
Parser.prototype.itf8 = function(e, t) {
  return this.setNextParser("itf8", e, t);
};
Parser.prototype.itf8 = function(e, t) {
  return this.setNextParser("itf8", e, t);
};
Parser.prototype.generateItf8 = function(e) {
  var t = e.generateVariable(this.varName), r = e.generateTmpVariable();
  e.pushCode(`
    var `.concat(
    r,
    ` = buffer[offset];
    if (`
  ).concat(
    r,
    ` < 0x80) {
      `
  ).concat(
    t,
    " = "
  ).concat(r, `;
      offset += 1;
    } else if (`).concat(
    r,
    ` < 0xc0) {
      `
  ).concat(
    t,
    " = (("
  ).concat(r, `<<8) | buffer[offset+1]) & 0x3fff;
      offset += 2;
    } else if (`).concat(
    r,
    ` < 0xe0) {
      `
  ).concat(
    t,
    " = (("
  ).concat(r, `<<16) | (buffer[offset+1]<< 8) |  buffer[offset+2]) & 0x1fffff;
      offset += 3;
    } else if (`).concat(
    r,
    ` < 0xf0) {
      `
  ).concat(
    t,
    " = (("
  ).concat(r, `<<24) | (buffer[offset+1]<<16) | (buffer[offset+2]<<8) | buffer[offset+3]) & 0x0fffffff;
      offset += 4
    } else {
      `).concat(
    t,
    " = (("
  ).concat(r, ` & 0x0f)<<28) | (buffer[offset+1]<<20) | (buffer[offset+2]<<12) | (buffer[offset+3]<<4) | (buffer[offset+4] & 0x0f);
      // x=((0xff & 0x0f)<<28) | (0xff<<20) | (0xff<<12) | (0xff<<4) | (0x0f & 0x0f);
      // TODO *val_p = uv < 0x80000000UL ? uv : -((int32_t) (0xffffffffUL - uv)) - 1;
      offset += 5
    }
  `));
};
Parser.prototype.ltf8 = function(e, t) {
  return this.setNextParser("ltf8", e, t);
};
Parser.prototype.generateLtf8 = function(e) {
  var t = e.generateVariable(this.varName), r = e.generateTmpVariable();
  e.pushCode(`
  var `.concat(
    r,
    ` = buffer[offset];
  if (`
  ).concat(
    r,
    ` < 0x80) {
    `
  ).concat(
    t,
    " = "
  ).concat(r, `;
    offset += 1;
  } else if (`).concat(
    r,
    ` < 0xc0) {
    `
  ).concat(
    t,
    ` = ((buffer[offset]<<8) | buffer[offset+1]) & 0x3fff;
    offset += 2;
  } else if (`
  ).concat(
    r,
    ` < 0xe0) {
    `
  ).concat(
    t,
    ` = ((buffer[offset]<<16) | (buffer[offset+1]<<8) | buffer[offset+2]) & 0x1fffff;
    `
  ).concat(
    t,
    " = ((("
  ).concat(r, ` & 63) << 16) | buffer.readUInt16LE(offset + 1));
    offset += 3;
  } else if (`).concat(
    r,
    ` < 0xf0) {
    `
  ).concat(
    t,
    ` = ((buffer[offset]<<24) | (buffer[offset+1]<<16) | (buffer[offset+2]<<8) | buffer[offset+3]) & 0x0fffffff;
    offset += 4;
  } else if (`
  ).concat(
    r,
    ` < 0xf8) {
    `
  ).concat(
    t,
    ` = (((buffer[offset] & 15) * Math.pow(2,32))) +
      (buffer[offset+1]<<24) | (buffer[offset+2]<<16 | buffer[offset+3]<<8 | buffer[offset+4])
    // TODO *val_p = uv < 0x80000000UL ? uv : -((int32_t) (0xffffffffUL - uv)) - 1;
    offset += 5;
  } else if (`
  ).concat(
    r,
    ` < 0xfc) {
    `
  ).concat(
    t,
    ` = ((((buffer[offset] & 7) << 8) | buffer[offset+1] )) * Math.pow(2,32) +
      (buffer[offset+2]<<24) | (buffer[offset+3]<<16 | buffer[offset+4]<<8 | buffer[offset+5])
    offset += 6;
  } else if (`
  ).concat(
    r,
    ` < 0xfe) {
    `
  ).concat(
    t,
    ` = ((((buffer[offset] & 3) << 16) | buffer[offset+1]<<8 | buffer[offset+2])) * Math.pow(2,32) +
      (buffer[offset+3]<<24) | (buffer[offset+4]<<16 | buffer[offset+5]<<8 | buffer[offset+6])
    offset += 7;
  } else if (`
  ).concat(
    r,
    ` < 0xff) {
    `
  ).concat(
    t,
    ` = Long.fromBytesBE(buffer.slice(offset+1,offset+8));
    if (`
  ).concat(
    t,
    ".greaterThan(Number.MAX_SAFE_INTEGER) || "
  ).concat(t, `.lessThan(Number.MIN_SAFE_INTEGER))
      throw new Error('integer overflow')
    `).concat(
    t,
    " = "
  ).concat(t, `.toNumber()
    offset += 8;
  } else {
    `).concat(
    t,
    ` = Long.fromBytesBE(buffer.slice(offset+1,offset+9));
    if (`
  ).concat(
    t,
    ".greaterThan(Number.MAX_SAFE_INTEGER) || "
  ).concat(t, `.lessThan(Number.MIN_SAFE_INTEGER))
      throw new Error('integer overflow')
    `).concat(
    t,
    " = "
  ).concat(t, `.toNumber()
    offset += 9;
  }
  `));
};
var Parser_1 = Parser;
export {
  Parser_1 as P
};
