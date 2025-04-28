import { a1 as ie, aJ as q, b3 as ge, Q as _t } from "./index-8ee3baa4.js";
var ar = { exports: {} }, se = typeof Reflect == "object" ? Reflect : null, hr = se && typeof se.apply == "function" ? se.apply : function(n, i, d) {
  return Function.prototype.apply.call(n, i, d);
}, he;
se && typeof se.ownKeys == "function" ? he = se.ownKeys : Object.getOwnPropertySymbols ? he = function(n) {
  return Object.getOwnPropertyNames(n).concat(Object.getOwnPropertySymbols(n));
} : he = function(n) {
  return Object.getOwnPropertyNames(n);
};
function Et(e) {
  console && console.warn && console.warn(e);
}
var tt = Number.isNaN || function(n) {
  return n !== n;
};
function z() {
  z.init.call(this);
}
ar.exports = z;
ar.exports.once = Tt;
z.EventEmitter = z;
z.prototype._events = void 0;
z.prototype._eventsCount = 0;
z.prototype._maxListeners = void 0;
var gr = 10;
function ve(e) {
  if (typeof e != "function")
    throw new TypeError('The "listener" argument must be of type Function. Received type ' + typeof e);
}
Object.defineProperty(z, "defaultMaxListeners", {
  enumerable: !0,
  get: function() {
    return gr;
  },
  set: function(e) {
    if (typeof e != "number" || e < 0 || tt(e))
      throw new RangeError('The value of "defaultMaxListeners" is out of range. It must be a non-negative number. Received ' + e + ".");
    gr = e;
  }
});
z.init = function() {
  (this._events === void 0 || this._events === Object.getPrototypeOf(this)._events) && (this._events = /* @__PURE__ */ Object.create(null), this._eventsCount = 0), this._maxListeners = this._maxListeners || void 0;
};
z.prototype.setMaxListeners = function(n) {
  if (typeof n != "number" || n < 0 || tt(n))
    throw new RangeError('The value of "n" is out of range. It must be a non-negative number. Received ' + n + ".");
  return this._maxListeners = n, this;
};
function nt(e) {
  return e._maxListeners === void 0 ? z.defaultMaxListeners : e._maxListeners;
}
z.prototype.getMaxListeners = function() {
  return nt(this);
};
z.prototype.emit = function(n) {
  for (var i = [], d = 1; d < arguments.length; d++)
    i.push(arguments[d]);
  var y = n === "error", s = this._events;
  if (s !== void 0)
    y = y && s.error === void 0;
  else if (!y)
    return !1;
  if (y) {
    var f;
    if (i.length > 0 && (f = i[0]), f instanceof Error)
      throw f;
    var l = new Error("Unhandled error." + (f ? " (" + f.message + ")" : ""));
    throw l.context = f, l;
  }
  var o = s[n];
  if (o === void 0)
    return !1;
  if (typeof o == "function")
    hr(o, this, i);
  else
    for (var h = o.length, v = ut(o, h), d = 0; d < h; ++d)
      hr(v[d], this, i);
  return !0;
};
function it(e, n, i, d) {
  var y, s, f;
  if (ve(i), s = e._events, s === void 0 ? (s = e._events = /* @__PURE__ */ Object.create(null), e._eventsCount = 0) : (s.newListener !== void 0 && (e.emit(
    "newListener",
    n,
    i.listener ? i.listener : i
  ), s = e._events), f = s[n]), f === void 0)
    f = s[n] = i, ++e._eventsCount;
  else if (typeof f == "function" ? f = s[n] = d ? [i, f] : [f, i] : d ? f.unshift(i) : f.push(i), y = nt(e), y > 0 && f.length > y && !f.warned) {
    f.warned = !0;
    var l = new Error("Possible EventEmitter memory leak detected. " + f.length + " " + String(n) + " listeners added. Use emitter.setMaxListeners() to increase limit");
    l.name = "MaxListenersExceededWarning", l.emitter = e, l.type = n, l.count = f.length, Et(l);
  }
  return e;
}
z.prototype.addListener = function(n, i) {
  return it(this, n, i, !1);
};
z.prototype.on = z.prototype.addListener;
z.prototype.prependListener = function(n, i) {
  return it(this, n, i, !0);
};
function At() {
  if (!this.fired)
    return this.target.removeListener(this.type, this.wrapFn), this.fired = !0, arguments.length === 0 ? this.listener.call(this.target) : this.listener.apply(this.target, arguments);
}
function ot(e, n, i) {
  var d = { fired: !1, wrapFn: void 0, target: e, type: n, listener: i }, y = At.bind(d);
  return y.listener = i, d.wrapFn = y, y;
}
z.prototype.once = function(n, i) {
  return ve(i), this.on(n, ot(this, n, i)), this;
};
z.prototype.prependOnceListener = function(n, i) {
  return ve(i), this.prependListener(n, ot(this, n, i)), this;
};
z.prototype.removeListener = function(n, i) {
  var d, y, s, f, l;
  if (ve(i), y = this._events, y === void 0)
    return this;
  if (d = y[n], d === void 0)
    return this;
  if (d === i || d.listener === i)
    --this._eventsCount === 0 ? this._events = /* @__PURE__ */ Object.create(null) : (delete y[n], y.removeListener && this.emit("removeListener", n, d.listener || i));
  else if (typeof d != "function") {
    for (s = -1, f = d.length - 1; f >= 0; f--)
      if (d[f] === i || d[f].listener === i) {
        l = d[f].listener, s = f;
        break;
      }
    if (s < 0)
      return this;
    s === 0 ? d.shift() : Rt(d, s), d.length === 1 && (y[n] = d[0]), y.removeListener !== void 0 && this.emit("removeListener", n, l || i);
  }
  return this;
};
z.prototype.off = z.prototype.removeListener;
z.prototype.removeAllListeners = function(n) {
  var i, d, y;
  if (d = this._events, d === void 0)
    return this;
  if (d.removeListener === void 0)
    return arguments.length === 0 ? (this._events = /* @__PURE__ */ Object.create(null), this._eventsCount = 0) : d[n] !== void 0 && (--this._eventsCount === 0 ? this._events = /* @__PURE__ */ Object.create(null) : delete d[n]), this;
  if (arguments.length === 0) {
    var s = Object.keys(d), f;
    for (y = 0; y < s.length; ++y)
      f = s[y], f !== "removeListener" && this.removeAllListeners(f);
    return this.removeAllListeners("removeListener"), this._events = /* @__PURE__ */ Object.create(null), this._eventsCount = 0, this;
  }
  if (i = d[n], typeof i == "function")
    this.removeListener(n, i);
  else if (i !== void 0)
    for (y = i.length - 1; y >= 0; y--)
      this.removeListener(n, i[y]);
  return this;
};
function at(e, n, i) {
  var d = e._events;
  if (d === void 0)
    return [];
  var y = d[n];
  return y === void 0 ? [] : typeof y == "function" ? i ? [y.listener || y] : [y] : i ? Ot(y) : ut(y, y.length);
}
z.prototype.listeners = function(n) {
  return at(this, n, !0);
};
z.prototype.rawListeners = function(n) {
  return at(this, n, !1);
};
z.listenerCount = function(e, n) {
  return typeof e.listenerCount == "function" ? e.listenerCount(n) : ft.call(e, n);
};
z.prototype.listenerCount = ft;
function ft(e) {
  var n = this._events;
  if (n !== void 0) {
    var i = n[e];
    if (typeof i == "function")
      return 1;
    if (i !== void 0)
      return i.length;
  }
  return 0;
}
z.prototype.eventNames = function() {
  return this._eventsCount > 0 ? he(this._events) : [];
};
function ut(e, n) {
  for (var i = new Array(n), d = 0; d < n; ++d)
    i[d] = e[d];
  return i;
}
function Rt(e, n) {
  for (; n + 1 < e.length; n++)
    e[n] = e[n + 1];
  e.pop();
}
function Ot(e) {
  for (var n = new Array(e.length), i = 0; i < n.length; ++i)
    n[i] = e[i].listener || e[i];
  return n;
}
function Tt(e, n) {
  return new Promise(function(i, d) {
    function y(f) {
      e.removeListener(n, s), d(f);
    }
    function s() {
      typeof e.removeListener == "function" && e.removeListener("error", y), i([].slice.call(arguments));
    }
    st(e, n, s, { once: !0 }), n !== "error" && Pt(e, y, { once: !0 });
  });
}
function Pt(e, n, i) {
  typeof e.on == "function" && st(e, "error", n, i);
}
function st(e, n, i, d) {
  if (typeof e.on == "function")
    d.once ? e.once(n, i) : e.on(n, i);
  else if (typeof e.addEventListener == "function")
    e.addEventListener(n, function y(s) {
      d.once && e.removeEventListener(n, y), i(s);
    });
  else
    throw new TypeError('The "emitter" argument must be of type EventEmitter. Received type ' + typeof e);
}
var fr = ar.exports, or = { exports: {} };
typeof Object.create == "function" ? or.exports = function(n, i) {
  i && (n.super_ = i, n.prototype = Object.create(i.prototype, {
    constructor: {
      value: n,
      enumerable: !1,
      writable: !0,
      configurable: !0
    }
  }));
} : or.exports = function(n, i) {
  if (i) {
    n.super_ = i;
    var d = function() {
    };
    d.prototype = i.prototype, n.prototype = new d(), n.prototype.constructor = n;
  }
};
var ae = or.exports, _e, vr;
function lt() {
  return vr || (vr = 1, _e = fr.EventEmitter), _e;
}
var Ee = {}, Ae = {}, Re, br;
function ct() {
  return br || (br = 1, Re = function() {
    if (typeof Symbol != "function" || typeof Object.getOwnPropertySymbols != "function")
      return !1;
    if (typeof Symbol.iterator == "symbol")
      return !0;
    var n = {}, i = Symbol("test"), d = Object(i);
    if (typeof i == "string" || Object.prototype.toString.call(i) !== "[object Symbol]" || Object.prototype.toString.call(d) !== "[object Symbol]")
      return !1;
    var y = 42;
    n[i] = y;
    for (i in n)
      return !1;
    if (typeof Object.keys == "function" && Object.keys(n).length !== 0 || typeof Object.getOwnPropertyNames == "function" && Object.getOwnPropertyNames(n).length !== 0)
      return !1;
    var s = Object.getOwnPropertySymbols(n);
    if (s.length !== 1 || s[0] !== i || !Object.prototype.propertyIsEnumerable.call(n, i))
      return !1;
    if (typeof Object.getOwnPropertyDescriptor == "function") {
      var f = Object.getOwnPropertyDescriptor(n, i);
      if (f.value !== y || f.enumerable !== !0)
        return !1;
    }
    return !0;
  }), Re;
}
var Oe, mr;
function be() {
  if (mr)
    return Oe;
  mr = 1;
  var e = ct();
  return Oe = function() {
    return e() && !!Symbol.toStringTag;
  }, Oe;
}
var Te, wr;
function jt() {
  if (wr)
    return Te;
  wr = 1;
  var e = typeof Symbol < "u" && Symbol, n = ct();
  return Te = function() {
    return typeof e != "function" || typeof Symbol != "function" || typeof e("foo") != "symbol" || typeof Symbol("bar") != "symbol" ? !1 : n();
  }, Te;
}
var Pe, Sr;
function Mt() {
  if (Sr)
    return Pe;
  Sr = 1;
  var e = {
    foo: {}
  }, n = Object;
  return Pe = function() {
    return { __proto__: e }.foo === e.foo && !({ __proto__: null } instanceof n);
  }, Pe;
}
var je, _r;
function Lt() {
  if (_r)
    return je;
  _r = 1;
  var e = "Function.prototype.bind called on incompatible ", n = Array.prototype.slice, i = Object.prototype.toString, d = "[object Function]";
  return je = function(s) {
    var f = this;
    if (typeof f != "function" || i.call(f) !== d)
      throw new TypeError(e + f);
    for (var l = n.call(arguments, 1), o, h = function() {
      if (this instanceof o) {
        var _ = f.apply(
          this,
          l.concat(n.call(arguments))
        );
        return Object(_) === _ ? _ : this;
      } else
        return f.apply(
          s,
          l.concat(n.call(arguments))
        );
    }, v = Math.max(0, f.length - l.length), S = [], O = 0; O < v; O++)
      S.push("$" + O);
    if (o = Function("binder", "return function (" + S.join(",") + "){ return binder.apply(this,arguments); }")(h), f.prototype) {
      var L = function() {
      };
      L.prototype = f.prototype, o.prototype = new L(), L.prototype = null;
    }
    return o;
  }, je;
}
var Me, Er;
function ur() {
  if (Er)
    return Me;
  Er = 1;
  var e = Lt();
  return Me = Function.prototype.bind || e, Me;
}
var Le, Ar;
function It() {
  if (Ar)
    return Le;
  Ar = 1;
  var e = ur();
  return Le = e.call(Function.call, Object.prototype.hasOwnProperty), Le;
}
var Ie, Rr;
function sr() {
  if (Rr)
    return Ie;
  Rr = 1;
  var e, n = SyntaxError, i = Function, d = TypeError, y = function(K) {
    try {
      return i('"use strict"; return (' + K + ").constructor;")();
    } catch {
    }
  }, s = Object.getOwnPropertyDescriptor;
  if (s)
    try {
      s({}, "");
    } catch {
      s = null;
    }
  var f = function() {
    throw new d();
  }, l = s ? function() {
    try {
      return arguments.callee, f;
    } catch {
      try {
        return s(arguments, "callee").get;
      } catch {
        return f;
      }
    }
  }() : f, o = jt()(), h = Mt()(), v = Object.getPrototypeOf || (h ? function(K) {
    return K.__proto__;
  } : null), S = {}, O = typeof Uint8Array > "u" || !v ? e : v(Uint8Array), L = {
    "%AggregateError%": typeof AggregateError > "u" ? e : AggregateError,
    "%Array%": Array,
    "%ArrayBuffer%": typeof ArrayBuffer > "u" ? e : ArrayBuffer,
    "%ArrayIteratorPrototype%": o && v ? v([][Symbol.iterator]()) : e,
    "%AsyncFromSyncIteratorPrototype%": e,
    "%AsyncFunction%": S,
    "%AsyncGenerator%": S,
    "%AsyncGeneratorFunction%": S,
    "%AsyncIteratorPrototype%": S,
    "%Atomics%": typeof Atomics > "u" ? e : Atomics,
    "%BigInt%": typeof BigInt > "u" ? e : BigInt,
    "%BigInt64Array%": typeof BigInt64Array > "u" ? e : BigInt64Array,
    "%BigUint64Array%": typeof BigUint64Array > "u" ? e : BigUint64Array,
    "%Boolean%": Boolean,
    "%DataView%": typeof DataView > "u" ? e : DataView,
    "%Date%": Date,
    "%decodeURI%": decodeURI,
    "%decodeURIComponent%": decodeURIComponent,
    "%encodeURI%": encodeURI,
    "%encodeURIComponent%": encodeURIComponent,
    "%Error%": Error,
    "%eval%": eval,
    // eslint-disable-line no-eval
    "%EvalError%": EvalError,
    "%Float32Array%": typeof Float32Array > "u" ? e : Float32Array,
    "%Float64Array%": typeof Float64Array > "u" ? e : Float64Array,
    "%FinalizationRegistry%": typeof FinalizationRegistry > "u" ? e : FinalizationRegistry,
    "%Function%": i,
    "%GeneratorFunction%": S,
    "%Int8Array%": typeof Int8Array > "u" ? e : Int8Array,
    "%Int16Array%": typeof Int16Array > "u" ? e : Int16Array,
    "%Int32Array%": typeof Int32Array > "u" ? e : Int32Array,
    "%isFinite%": isFinite,
    "%isNaN%": isNaN,
    "%IteratorPrototype%": o && v ? v(v([][Symbol.iterator]())) : e,
    "%JSON%": typeof JSON == "object" ? JSON : e,
    "%Map%": typeof Map > "u" ? e : Map,
    "%MapIteratorPrototype%": typeof Map > "u" || !o || !v ? e : v((/* @__PURE__ */ new Map())[Symbol.iterator]()),
    "%Math%": Math,
    "%Number%": Number,
    "%Object%": Object,
    "%parseFloat%": parseFloat,
    "%parseInt%": parseInt,
    "%Promise%": typeof Promise > "u" ? e : Promise,
    "%Proxy%": typeof Proxy > "u" ? e : Proxy,
    "%RangeError%": RangeError,
    "%ReferenceError%": ReferenceError,
    "%Reflect%": typeof Reflect > "u" ? e : Reflect,
    "%RegExp%": RegExp,
    "%Set%": typeof Set > "u" ? e : Set,
    "%SetIteratorPrototype%": typeof Set > "u" || !o || !v ? e : v((/* @__PURE__ */ new Set())[Symbol.iterator]()),
    "%SharedArrayBuffer%": typeof SharedArrayBuffer > "u" ? e : SharedArrayBuffer,
    "%String%": String,
    "%StringIteratorPrototype%": o && v ? v(""[Symbol.iterator]()) : e,
    "%Symbol%": o ? Symbol : e,
    "%SyntaxError%": n,
    "%ThrowTypeError%": l,
    "%TypedArray%": O,
    "%TypeError%": d,
    "%Uint8Array%": typeof Uint8Array > "u" ? e : Uint8Array,
    "%Uint8ClampedArray%": typeof Uint8ClampedArray > "u" ? e : Uint8ClampedArray,
    "%Uint16Array%": typeof Uint16Array > "u" ? e : Uint16Array,
    "%Uint32Array%": typeof Uint32Array > "u" ? e : Uint32Array,
    "%URIError%": URIError,
    "%WeakMap%": typeof WeakMap > "u" ? e : WeakMap,
    "%WeakRef%": typeof WeakRef > "u" ? e : WeakRef,
    "%WeakSet%": typeof WeakSet > "u" ? e : WeakSet
  };
  if (v)
    try {
      null.error;
    } catch (K) {
      var _ = v(v(K));
      L["%Error.prototype%"] = _;
    }
  var j = function K(W) {
    var $;
    if (W === "%AsyncFunction%")
      $ = y("async function () {}");
    else if (W === "%GeneratorFunction%")
      $ = y("function* () {}");
    else if (W === "%AsyncGeneratorFunction%")
      $ = y("async function* () {}");
    else if (W === "%AsyncGenerator%") {
      var G = K("%AsyncGeneratorFunction%");
      G && ($ = G.prototype);
    } else if (W === "%AsyncIteratorPrototype%") {
      var Y = K("%AsyncGenerator%");
      Y && v && ($ = v(Y.prototype));
    }
    return L[W] = $, $;
  }, E = {
    "%ArrayBufferPrototype%": ["ArrayBuffer", "prototype"],
    "%ArrayPrototype%": ["Array", "prototype"],
    "%ArrayProto_entries%": ["Array", "prototype", "entries"],
    "%ArrayProto_forEach%": ["Array", "prototype", "forEach"],
    "%ArrayProto_keys%": ["Array", "prototype", "keys"],
    "%ArrayProto_values%": ["Array", "prototype", "values"],
    "%AsyncFunctionPrototype%": ["AsyncFunction", "prototype"],
    "%AsyncGenerator%": ["AsyncGeneratorFunction", "prototype"],
    "%AsyncGeneratorPrototype%": ["AsyncGeneratorFunction", "prototype", "prototype"],
    "%BooleanPrototype%": ["Boolean", "prototype"],
    "%DataViewPrototype%": ["DataView", "prototype"],
    "%DatePrototype%": ["Date", "prototype"],
    "%ErrorPrototype%": ["Error", "prototype"],
    "%EvalErrorPrototype%": ["EvalError", "prototype"],
    "%Float32ArrayPrototype%": ["Float32Array", "prototype"],
    "%Float64ArrayPrototype%": ["Float64Array", "prototype"],
    "%FunctionPrototype%": ["Function", "prototype"],
    "%Generator%": ["GeneratorFunction", "prototype"],
    "%GeneratorPrototype%": ["GeneratorFunction", "prototype", "prototype"],
    "%Int8ArrayPrototype%": ["Int8Array", "prototype"],
    "%Int16ArrayPrototype%": ["Int16Array", "prototype"],
    "%Int32ArrayPrototype%": ["Int32Array", "prototype"],
    "%JSONParse%": ["JSON", "parse"],
    "%JSONStringify%": ["JSON", "stringify"],
    "%MapPrototype%": ["Map", "prototype"],
    "%NumberPrototype%": ["Number", "prototype"],
    "%ObjectPrototype%": ["Object", "prototype"],
    "%ObjProto_toString%": ["Object", "prototype", "toString"],
    "%ObjProto_valueOf%": ["Object", "prototype", "valueOf"],
    "%PromisePrototype%": ["Promise", "prototype"],
    "%PromiseProto_then%": ["Promise", "prototype", "then"],
    "%Promise_all%": ["Promise", "all"],
    "%Promise_reject%": ["Promise", "reject"],
    "%Promise_resolve%": ["Promise", "resolve"],
    "%RangeErrorPrototype%": ["RangeError", "prototype"],
    "%ReferenceErrorPrototype%": ["ReferenceError", "prototype"],
    "%RegExpPrototype%": ["RegExp", "prototype"],
    "%SetPrototype%": ["Set", "prototype"],
    "%SharedArrayBufferPrototype%": ["SharedArrayBuffer", "prototype"],
    "%StringPrototype%": ["String", "prototype"],
    "%SymbolPrototype%": ["Symbol", "prototype"],
    "%SyntaxErrorPrototype%": ["SyntaxError", "prototype"],
    "%TypedArrayPrototype%": ["TypedArray", "prototype"],
    "%TypeErrorPrototype%": ["TypeError", "prototype"],
    "%Uint8ArrayPrototype%": ["Uint8Array", "prototype"],
    "%Uint8ClampedArrayPrototype%": ["Uint8ClampedArray", "prototype"],
    "%Uint16ArrayPrototype%": ["Uint16Array", "prototype"],
    "%Uint32ArrayPrototype%": ["Uint32Array", "prototype"],
    "%URIErrorPrototype%": ["URIError", "prototype"],
    "%WeakMapPrototype%": ["WeakMap", "prototype"],
    "%WeakSetPrototype%": ["WeakSet", "prototype"]
  }, a = ur(), g = It(), b = a.call(Function.call, Array.prototype.concat), P = a.call(Function.apply, Array.prototype.splice), I = a.call(Function.call, String.prototype.replace), D = a.call(Function.call, String.prototype.slice), x = a.call(Function.call, RegExp.prototype.exec), U = /[^%.[\]]+|\[(?:(-?\d+(?:\.\d+)?)|(["'])((?:(?!\2)[^\\]|\\.)*?)\2)\]|(?=(?:\.|\[\])(?:\.|\[\]|%$))/g, C = /\\(\\)?/g, Q = function(W) {
    var $ = D(W, 0, 1), G = D(W, -1);
    if ($ === "%" && G !== "%")
      throw new n("invalid intrinsic syntax, expected closing `%`");
    if (G === "%" && $ !== "%")
      throw new n("invalid intrinsic syntax, expected opening `%`");
    var Y = [];
    return I(W, U, function(Z, X, H, ee) {
      Y[Y.length] = H ? I(ee, C, "$1") : X || Z;
    }), Y;
  }, re = function(W, $) {
    var G = W, Y;
    if (g(E, G) && (Y = E[G], G = "%" + Y[0] + "%"), g(L, G)) {
      var Z = L[G];
      if (Z === S && (Z = j(G)), typeof Z > "u" && !$)
        throw new d("intrinsic " + W + " exists, but is not available. Please file an issue!");
      return {
        alias: Y,
        name: G,
        value: Z
      };
    }
    throw new n("intrinsic " + W + " does not exist!");
  };
  return Ie = function(W, $) {
    if (typeof W != "string" || W.length === 0)
      throw new d("intrinsic name must be a non-empty string");
    if (arguments.length > 1 && typeof $ != "boolean")
      throw new d('"allowMissing" argument must be a boolean');
    if (x(/^%?[^%]*%?$/, W) === null)
      throw new n("`%` may not be present anywhere but at the beginning and end of the intrinsic name");
    var G = Q(W), Y = G.length > 0 ? G[0] : "", Z = re("%" + Y + "%", $), X = Z.name, H = Z.value, ee = !1, te = Z.alias;
    te && (Y = te[0], P(G, b([0, 1], te)));
    for (var p = 1, m = !0; p < G.length; p += 1) {
      var R = G[p], F = D(R, 0, 1), c = D(R, -1);
      if ((F === '"' || F === "'" || F === "`" || c === '"' || c === "'" || c === "`") && F !== c)
        throw new n("property names with quotes must have matching quotes");
      if ((R === "constructor" || !m) && (ee = !0), Y += "." + R, X = "%" + Y + "%", g(L, X))
        H = L[X];
      else if (H != null) {
        if (!(R in H)) {
          if (!$)
            throw new d("base intrinsic for " + W + " exists, but the property is not available.");
          return;
        }
        if (s && p + 1 >= G.length) {
          var u = s(H, R);
          m = !!u, m && "get" in u && !("originalValue" in u.get) ? H = u.get : H = H[R];
        } else
          m = g(H, R), H = H[R];
        m && !ee && (L[X] = H);
      }
    }
    return H;
  }, Ie;
}
var Ne = { exports: {} }, Or;
function Nt() {
  return Or || (Or = 1, function(e) {
    var n = ur(), i = sr(), d = i("%Function.prototype.apply%"), y = i("%Function.prototype.call%"), s = i("%Reflect.apply%", !0) || n.call(y, d), f = i("%Object.getOwnPropertyDescriptor%", !0), l = i("%Object.defineProperty%", !0), o = i("%Math.max%");
    if (l)
      try {
        l({}, "a", { value: 1 });
      } catch {
        l = null;
      }
    e.exports = function(S) {
      var O = s(n, y, arguments);
      if (f && l) {
        var L = f(O, "length");
        L.configurable && l(
          O,
          "length",
          { value: 1 + o(0, S.length - (arguments.length - 1)) }
        );
      }
      return O;
    };
    var h = function() {
      return s(n, d, arguments);
    };
    l ? l(e.exports, "apply", { value: h }) : e.exports.apply = h;
  }(Ne)), Ne.exports;
}
var Be, Tr;
function lr() {
  if (Tr)
    return Be;
  Tr = 1;
  var e = sr(), n = Nt(), i = n(e("String.prototype.indexOf"));
  return Be = function(y, s) {
    var f = e(y, !!s);
    return typeof f == "function" && i(y, ".prototype.") > -1 ? n(f) : f;
  }, Be;
}
var Ce, Pr;
function Bt() {
  if (Pr)
    return Ce;
  Pr = 1;
  var e = be()(), n = lr(), i = n("Object.prototype.toString"), d = function(l) {
    return e && l && typeof l == "object" && Symbol.toStringTag in l ? !1 : i(l) === "[object Arguments]";
  }, y = function(l) {
    return d(l) ? !0 : l !== null && typeof l == "object" && typeof l.length == "number" && l.length >= 0 && i(l) !== "[object Array]" && i(l.callee) === "[object Function]";
  }, s = function() {
    return d(arguments);
  }();
  return d.isLegacyArguments = y, Ce = s ? d : y, Ce;
}
var De, jr;
function Ct() {
  if (jr)
    return De;
  jr = 1;
  var e = Object.prototype.toString, n = Function.prototype.toString, i = /^\s*(?:function)?\*/, d = be()(), y = Object.getPrototypeOf, s = function() {
    if (!d)
      return !1;
    try {
      return Function("return function*() {}")();
    } catch {
    }
  }, f;
  return De = function(o) {
    if (typeof o != "function")
      return !1;
    if (i.test(n.call(o)))
      return !0;
    if (!d) {
      var h = e.call(o);
      return h === "[object GeneratorFunction]";
    }
    if (!y)
      return !1;
    if (typeof f > "u") {
      var v = s();
      f = v ? y(v) : !1;
    }
    return y(o) === f;
  }, De;
}
var Fe, Mr;
function Dt() {
  if (Mr)
    return Fe;
  Mr = 1;
  var e = Function.prototype.toString, n = typeof Reflect == "object" && Reflect !== null && Reflect.apply, i, d;
  if (typeof n == "function" && typeof Object.defineProperty == "function")
    try {
      i = Object.defineProperty({}, "length", {
        get: function() {
          throw d;
        }
      }), d = {}, n(function() {
        throw 42;
      }, null, i);
    } catch (g) {
      g !== d && (n = null);
    }
  else
    n = null;
  var y = /^\s*class\b/, s = function(b) {
    try {
      var P = e.call(b);
      return y.test(P);
    } catch {
      return !1;
    }
  }, f = function(b) {
    try {
      return s(b) ? !1 : (e.call(b), !0);
    } catch {
      return !1;
    }
  }, l = Object.prototype.toString, o = "[object Object]", h = "[object Function]", v = "[object GeneratorFunction]", S = "[object HTMLAllCollection]", O = "[object HTML document.all class]", L = "[object HTMLCollection]", _ = typeof Symbol == "function" && !!Symbol.toStringTag, j = !(0 in [,]), E = function() {
    return !1;
  };
  if (typeof document == "object") {
    var a = document.all;
    l.call(a) === l.call(document.all) && (E = function(b) {
      if ((j || !b) && (typeof b > "u" || typeof b == "object"))
        try {
          var P = l.call(b);
          return (P === S || P === O || P === L || P === o) && b("") == null;
        } catch {
        }
      return !1;
    });
  }
  return Fe = n ? function(b) {
    if (E(b))
      return !0;
    if (!b || typeof b != "function" && typeof b != "object")
      return !1;
    try {
      n(b, null, i);
    } catch (P) {
      if (P !== d)
        return !1;
    }
    return !s(b) && f(b);
  } : function(b) {
    if (E(b))
      return !0;
    if (!b || typeof b != "function" && typeof b != "object")
      return !1;
    if (_)
      return f(b);
    if (s(b))
      return !1;
    var P = l.call(b);
    return P !== h && P !== v && !/^\[object HTML/.test(P) ? !1 : f(b);
  }, Fe;
}
var Ue, Lr;
function dt() {
  if (Lr)
    return Ue;
  Lr = 1;
  var e = Dt(), n = Object.prototype.toString, i = Object.prototype.hasOwnProperty, d = function(o, h, v) {
    for (var S = 0, O = o.length; S < O; S++)
      i.call(o, S) && (v == null ? h(o[S], S, o) : h.call(v, o[S], S, o));
  }, y = function(o, h, v) {
    for (var S = 0, O = o.length; S < O; S++)
      v == null ? h(o.charAt(S), S, o) : h.call(v, o.charAt(S), S, o);
  }, s = function(o, h, v) {
    for (var S in o)
      i.call(o, S) && (v == null ? h(o[S], S, o) : h.call(v, o[S], S, o));
  }, f = function(o, h, v) {
    if (!e(h))
      throw new TypeError("iterator must be a function");
    var S;
    arguments.length >= 3 && (S = v), n.call(o) === "[object Array]" ? d(o, h, S) : typeof o == "string" ? y(o, h, S) : s(o, h, S);
  };
  return Ue = f, Ue;
}
var ke, Ir;
function pt() {
  if (Ir)
    return ke;
  Ir = 1;
  var e = [
    "BigInt64Array",
    "BigUint64Array",
    "Float32Array",
    "Float64Array",
    "Int16Array",
    "Int32Array",
    "Int8Array",
    "Uint16Array",
    "Uint32Array",
    "Uint8Array",
    "Uint8ClampedArray"
  ], n = typeof globalThis > "u" ? ie : globalThis;
  return ke = function() {
    for (var d = [], y = 0; y < e.length; y++)
      typeof n[e[y]] == "function" && (d[d.length] = e[y]);
    return d;
  }, ke;
}
var qe, Nr;
function yt() {
  if (Nr)
    return qe;
  Nr = 1;
  var e = sr(), n = e("%Object.getOwnPropertyDescriptor%", !0);
  if (n)
    try {
      n([], "length");
    } catch {
      n = null;
    }
  return qe = n, qe;
}
var We, Br;
function ht() {
  if (Br)
    return We;
  Br = 1;
  var e = dt(), n = pt(), i = lr(), d = i("Object.prototype.toString"), y = be()(), s = yt(), f = typeof globalThis > "u" ? ie : globalThis, l = n(), o = i("Array.prototype.indexOf", !0) || function(_, j) {
    for (var E = 0; E < _.length; E += 1)
      if (_[E] === j)
        return E;
    return -1;
  }, h = i("String.prototype.slice"), v = {}, S = Object.getPrototypeOf;
  y && s && S && e(l, function(L) {
    var _ = new f[L]();
    if (Symbol.toStringTag in _) {
      var j = S(_), E = s(j, Symbol.toStringTag);
      if (!E) {
        var a = S(j);
        E = s(a, Symbol.toStringTag);
      }
      v[L] = E.get;
    }
  });
  var O = function(_) {
    var j = !1;
    return e(v, function(E, a) {
      if (!j)
        try {
          j = E.call(_) === a;
        } catch {
        }
    }), j;
  };
  return We = function(_) {
    if (!_ || typeof _ != "object")
      return !1;
    if (!y || !(Symbol.toStringTag in _)) {
      var j = h(d(_), 8, -1);
      return o(l, j) > -1;
    }
    return s ? O(_) : !1;
  }, We;
}
var xe, Cr;
function Ft() {
  if (Cr)
    return xe;
  Cr = 1;
  var e = dt(), n = pt(), i = lr(), d = yt(), y = i("Object.prototype.toString"), s = be()(), f = typeof globalThis > "u" ? ie : globalThis, l = n(), o = i("String.prototype.slice"), h = {}, v = Object.getPrototypeOf;
  s && d && v && e(l, function(L) {
    if (typeof f[L] == "function") {
      var _ = new f[L]();
      if (Symbol.toStringTag in _) {
        var j = v(_), E = d(j, Symbol.toStringTag);
        if (!E) {
          var a = v(j);
          E = d(a, Symbol.toStringTag);
        }
        h[L] = E.get;
      }
    }
  });
  var S = function(_) {
    var j = !1;
    return e(h, function(E, a) {
      if (!j)
        try {
          var g = E.call(_);
          g === a && (j = g);
        } catch {
        }
    }), j;
  }, O = ht();
  return xe = function(_) {
    return O(_) ? !s || !(Symbol.toStringTag in _) ? o(y(_), 8, -1) : S(_) : !1;
  }, xe;
}
var Dr;
function Ut() {
  return Dr || (Dr = 1, function(e) {
    var n = Bt(), i = Ct(), d = Ft(), y = ht();
    function s(A) {
      return A.call.bind(A);
    }
    var f = typeof BigInt < "u", l = typeof Symbol < "u", o = s(Object.prototype.toString), h = s(Number.prototype.valueOf), v = s(String.prototype.valueOf), S = s(Boolean.prototype.valueOf);
    if (f)
      var O = s(BigInt.prototype.valueOf);
    if (l)
      var L = s(Symbol.prototype.valueOf);
    function _(A, pe) {
      if (typeof A != "object")
        return !1;
      try {
        return pe(A), !0;
      } catch {
        return !1;
      }
    }
    e.isArgumentsObject = n, e.isGeneratorFunction = i, e.isTypedArray = y;
    function j(A) {
      return typeof Promise < "u" && A instanceof Promise || A !== null && typeof A == "object" && typeof A.then == "function" && typeof A.catch == "function";
    }
    e.isPromise = j;
    function E(A) {
      return typeof ArrayBuffer < "u" && ArrayBuffer.isView ? ArrayBuffer.isView(A) : y(A) || m(A);
    }
    e.isArrayBufferView = E;
    function a(A) {
      return d(A) === "Uint8Array";
    }
    e.isUint8Array = a;
    function g(A) {
      return d(A) === "Uint8ClampedArray";
    }
    e.isUint8ClampedArray = g;
    function b(A) {
      return d(A) === "Uint16Array";
    }
    e.isUint16Array = b;
    function P(A) {
      return d(A) === "Uint32Array";
    }
    e.isUint32Array = P;
    function I(A) {
      return d(A) === "Int8Array";
    }
    e.isInt8Array = I;
    function D(A) {
      return d(A) === "Int16Array";
    }
    e.isInt16Array = D;
    function x(A) {
      return d(A) === "Int32Array";
    }
    e.isInt32Array = x;
    function U(A) {
      return d(A) === "Float32Array";
    }
    e.isFloat32Array = U;
    function C(A) {
      return d(A) === "Float64Array";
    }
    e.isFloat64Array = C;
    function Q(A) {
      return d(A) === "BigInt64Array";
    }
    e.isBigInt64Array = Q;
    function re(A) {
      return d(A) === "BigUint64Array";
    }
    e.isBigUint64Array = re;
    function K(A) {
      return o(A) === "[object Map]";
    }
    K.working = typeof Map < "u" && K(/* @__PURE__ */ new Map());
    function W(A) {
      return typeof Map > "u" ? !1 : K.working ? K(A) : A instanceof Map;
    }
    e.isMap = W;
    function $(A) {
      return o(A) === "[object Set]";
    }
    $.working = typeof Set < "u" && $(/* @__PURE__ */ new Set());
    function G(A) {
      return typeof Set > "u" ? !1 : $.working ? $(A) : A instanceof Set;
    }
    e.isSet = G;
    function Y(A) {
      return o(A) === "[object WeakMap]";
    }
    Y.working = typeof WeakMap < "u" && Y(/* @__PURE__ */ new WeakMap());
    function Z(A) {
      return typeof WeakMap > "u" ? !1 : Y.working ? Y(A) : A instanceof WeakMap;
    }
    e.isWeakMap = Z;
    function X(A) {
      return o(A) === "[object WeakSet]";
    }
    X.working = typeof WeakSet < "u" && X(/* @__PURE__ */ new WeakSet());
    function H(A) {
      return X(A);
    }
    e.isWeakSet = H;
    function ee(A) {
      return o(A) === "[object ArrayBuffer]";
    }
    ee.working = typeof ArrayBuffer < "u" && ee(new ArrayBuffer());
    function te(A) {
      return typeof ArrayBuffer > "u" ? !1 : ee.working ? ee(A) : A instanceof ArrayBuffer;
    }
    e.isArrayBuffer = te;
    function p(A) {
      return o(A) === "[object DataView]";
    }
    p.working = typeof ArrayBuffer < "u" && typeof DataView < "u" && p(new DataView(new ArrayBuffer(1), 0, 1));
    function m(A) {
      return typeof DataView > "u" ? !1 : p.working ? p(A) : A instanceof DataView;
    }
    e.isDataView = m;
    var R = typeof SharedArrayBuffer < "u" ? SharedArrayBuffer : void 0;
    function F(A) {
      return o(A) === "[object SharedArrayBuffer]";
    }
    function c(A) {
      return typeof R > "u" ? !1 : (typeof F.working > "u" && (F.working = F(new R())), F.working ? F(A) : A instanceof R);
    }
    e.isSharedArrayBuffer = c;
    function u(A) {
      return o(A) === "[object AsyncFunction]";
    }
    e.isAsyncFunction = u;
    function w(A) {
      return o(A) === "[object Map Iterator]";
    }
    e.isMapIterator = w;
    function M(A) {
      return o(A) === "[object Set Iterator]";
    }
    e.isSetIterator = M;
    function B(A) {
      return o(A) === "[object Generator]";
    }
    e.isGeneratorObject = B;
    function r(A) {
      return o(A) === "[object WebAssembly.Module]";
    }
    e.isWebAssemblyCompiledModule = r;
    function t(A) {
      return _(A, h);
    }
    e.isNumberObject = t;
    function T(A) {
      return _(A, v);
    }
    e.isStringObject = T;
    function N(A) {
      return _(A, S);
    }
    e.isBooleanObject = N;
    function J(A) {
      return f && _(A, O);
    }
    e.isBigIntObject = J;
    function k(A) {
      return l && _(A, L);
    }
    e.isSymbolObject = k;
    function V(A) {
      return t(A) || T(A) || N(A) || J(A) || k(A);
    }
    e.isBoxedPrimitive = V;
    function oe(A) {
      return typeof Uint8Array < "u" && (te(A) || c(A));
    }
    e.isAnyArrayBuffer = oe, ["isProxy", "isExternal", "isModuleNamespaceObject"].forEach(function(A) {
      Object.defineProperty(e, A, {
        enumerable: !1,
        value: function() {
          throw new Error(A + " is not supported in userland");
        }
      });
    });
  }(Ae)), Ae;
}
var Ge, Fr;
function kt() {
  return Fr || (Fr = 1, Ge = function(n) {
    return n && typeof n == "object" && typeof n.copy == "function" && typeof n.fill == "function" && typeof n.readUInt8 == "function";
  }), Ge;
}
var Ur;
function gt() {
  return Ur || (Ur = 1, function(e) {
    var n = Object.getOwnPropertyDescriptors || function(m) {
      for (var R = Object.keys(m), F = {}, c = 0; c < R.length; c++)
        F[R[c]] = Object.getOwnPropertyDescriptor(m, R[c]);
      return F;
    }, i = /%[sdj%]/g;
    e.format = function(p) {
      if (!I(p)) {
        for (var m = [], R = 0; R < arguments.length; R++)
          m.push(f(arguments[R]));
        return m.join(" ");
      }
      for (var R = 1, F = arguments, c = F.length, u = String(p).replace(i, function(M) {
        if (M === "%%")
          return "%";
        if (R >= c)
          return M;
        switch (M) {
          case "%s":
            return String(F[R++]);
          case "%d":
            return Number(F[R++]);
          case "%j":
            try {
              return JSON.stringify(F[R++]);
            } catch {
              return "[Circular]";
            }
          default:
            return M;
        }
      }), w = F[R]; R < c; w = F[++R])
        g(w) || !C(w) ? u += " " + w : u += " " + f(w);
      return u;
    }, e.deprecate = function(p, m) {
      if (typeof q < "u" && q.noDeprecation === !0)
        return p;
      if (typeof q > "u")
        return function() {
          return e.deprecate(p, m).apply(this, arguments);
        };
      var R = !1;
      function F() {
        if (!R) {
          if (q.throwDeprecation)
            throw new Error(m);
          q.traceDeprecation ? console.trace(m) : console.error(m), R = !0;
        }
        return p.apply(this, arguments);
      }
      return F;
    };
    var d = {}, y = /^$/;
    if (q.env.NODE_DEBUG) {
      var s = q.env.NODE_DEBUG;
      s = s.replace(/[|\\{}()[\]^$+?.]/g, "\\$&").replace(/\*/g, ".*").replace(/,/g, "$|^").toUpperCase(), y = new RegExp("^" + s + "$", "i");
    }
    e.debuglog = function(p) {
      if (p = p.toUpperCase(), !d[p])
        if (y.test(p)) {
          var m = q.pid;
          d[p] = function() {
            var R = e.format.apply(e, arguments);
            console.error("%s %d: %s", p, m, R);
          };
        } else
          d[p] = function() {
          };
      return d[p];
    };
    function f(p, m) {
      var R = {
        seen: [],
        stylize: o
      };
      return arguments.length >= 3 && (R.depth = arguments[2]), arguments.length >= 4 && (R.colors = arguments[3]), a(m) ? R.showHidden = m : m && e._extend(R, m), x(R.showHidden) && (R.showHidden = !1), x(R.depth) && (R.depth = 2), x(R.colors) && (R.colors = !1), x(R.customInspect) && (R.customInspect = !0), R.colors && (R.stylize = l), v(R, p, R.depth);
    }
    e.inspect = f, f.colors = {
      bold: [1, 22],
      italic: [3, 23],
      underline: [4, 24],
      inverse: [7, 27],
      white: [37, 39],
      grey: [90, 39],
      black: [30, 39],
      blue: [34, 39],
      cyan: [36, 39],
      green: [32, 39],
      magenta: [35, 39],
      red: [31, 39],
      yellow: [33, 39]
    }, f.styles = {
      special: "cyan",
      number: "yellow",
      boolean: "yellow",
      undefined: "grey",
      null: "bold",
      string: "green",
      date: "magenta",
      // "name": intentionally not styling
      regexp: "red"
    };
    function l(p, m) {
      var R = f.styles[m];
      return R ? "\x1B[" + f.colors[R][0] + "m" + p + "\x1B[" + f.colors[R][1] + "m" : p;
    }
    function o(p, m) {
      return p;
    }
    function h(p) {
      var m = {};
      return p.forEach(function(R, F) {
        m[R] = !0;
      }), m;
    }
    function v(p, m, R) {
      if (p.customInspect && m && K(m.inspect) && // Filter out the util module, it's inspect function is special
      m.inspect !== e.inspect && // Also filter out any prototype objects using the circular check.
      !(m.constructor && m.constructor.prototype === m)) {
        var F = m.inspect(R, p);
        return I(F) || (F = v(p, F, R)), F;
      }
      var c = S(p, m);
      if (c)
        return c;
      var u = Object.keys(m), w = h(u);
      if (p.showHidden && (u = Object.getOwnPropertyNames(m)), re(m) && (u.indexOf("message") >= 0 || u.indexOf("description") >= 0))
        return O(m);
      if (u.length === 0) {
        if (K(m)) {
          var M = m.name ? ": " + m.name : "";
          return p.stylize("[Function" + M + "]", "special");
        }
        if (U(m))
          return p.stylize(RegExp.prototype.toString.call(m), "regexp");
        if (Q(m))
          return p.stylize(Date.prototype.toString.call(m), "date");
        if (re(m))
          return O(m);
      }
      var B = "", r = !1, t = ["{", "}"];
      if (E(m) && (r = !0, t = ["[", "]"]), K(m)) {
        var T = m.name ? ": " + m.name : "";
        B = " [Function" + T + "]";
      }
      if (U(m) && (B = " " + RegExp.prototype.toString.call(m)), Q(m) && (B = " " + Date.prototype.toUTCString.call(m)), re(m) && (B = " " + O(m)), u.length === 0 && (!r || m.length == 0))
        return t[0] + B + t[1];
      if (R < 0)
        return U(m) ? p.stylize(RegExp.prototype.toString.call(m), "regexp") : p.stylize("[Object]", "special");
      p.seen.push(m);
      var N;
      return r ? N = L(p, m, R, w, u) : N = u.map(function(J) {
        return _(p, m, R, w, J, r);
      }), p.seen.pop(), j(N, B, t);
    }
    function S(p, m) {
      if (x(m))
        return p.stylize("undefined", "undefined");
      if (I(m)) {
        var R = "'" + JSON.stringify(m).replace(/^"|"$/g, "").replace(/'/g, "\\'").replace(/\\"/g, '"') + "'";
        return p.stylize(R, "string");
      }
      if (P(m))
        return p.stylize("" + m, "number");
      if (a(m))
        return p.stylize("" + m, "boolean");
      if (g(m))
        return p.stylize("null", "null");
    }
    function O(p) {
      return "[" + Error.prototype.toString.call(p) + "]";
    }
    function L(p, m, R, F, c) {
      for (var u = [], w = 0, M = m.length; w < M; ++w)
        X(m, String(w)) ? u.push(_(
          p,
          m,
          R,
          F,
          String(w),
          !0
        )) : u.push("");
      return c.forEach(function(B) {
        B.match(/^\d+$/) || u.push(_(
          p,
          m,
          R,
          F,
          B,
          !0
        ));
      }), u;
    }
    function _(p, m, R, F, c, u) {
      var w, M, B;
      if (B = Object.getOwnPropertyDescriptor(m, c) || { value: m[c] }, B.get ? B.set ? M = p.stylize("[Getter/Setter]", "special") : M = p.stylize("[Getter]", "special") : B.set && (M = p.stylize("[Setter]", "special")), X(F, c) || (w = "[" + c + "]"), M || (p.seen.indexOf(B.value) < 0 ? (g(R) ? M = v(p, B.value, null) : M = v(p, B.value, R - 1), M.indexOf(`
`) > -1 && (u ? M = M.split(`
`).map(function(r) {
        return "  " + r;
      }).join(`
`).slice(2) : M = `
` + M.split(`
`).map(function(r) {
        return "   " + r;
      }).join(`
`))) : M = p.stylize("[Circular]", "special")), x(w)) {
        if (u && c.match(/^\d+$/))
          return M;
        w = JSON.stringify("" + c), w.match(/^"([a-zA-Z_][a-zA-Z_0-9]*)"$/) ? (w = w.slice(1, -1), w = p.stylize(w, "name")) : (w = w.replace(/'/g, "\\'").replace(/\\"/g, '"').replace(/(^"|"$)/g, "'"), w = p.stylize(w, "string"));
      }
      return w + ": " + M;
    }
    function j(p, m, R) {
      var F = p.reduce(function(c, u) {
        return u.indexOf(`
`) >= 0, c + u.replace(/\u001b\[\d\d?m/g, "").length + 1;
      }, 0);
      return F > 60 ? R[0] + (m === "" ? "" : m + `
 `) + " " + p.join(`,
  `) + " " + R[1] : R[0] + m + " " + p.join(", ") + " " + R[1];
    }
    e.types = Ut();
    function E(p) {
      return Array.isArray(p);
    }
    e.isArray = E;
    function a(p) {
      return typeof p == "boolean";
    }
    e.isBoolean = a;
    function g(p) {
      return p === null;
    }
    e.isNull = g;
    function b(p) {
      return p == null;
    }
    e.isNullOrUndefined = b;
    function P(p) {
      return typeof p == "number";
    }
    e.isNumber = P;
    function I(p) {
      return typeof p == "string";
    }
    e.isString = I;
    function D(p) {
      return typeof p == "symbol";
    }
    e.isSymbol = D;
    function x(p) {
      return p === void 0;
    }
    e.isUndefined = x;
    function U(p) {
      return C(p) && $(p) === "[object RegExp]";
    }
    e.isRegExp = U, e.types.isRegExp = U;
    function C(p) {
      return typeof p == "object" && p !== null;
    }
    e.isObject = C;
    function Q(p) {
      return C(p) && $(p) === "[object Date]";
    }
    e.isDate = Q, e.types.isDate = Q;
    function re(p) {
      return C(p) && ($(p) === "[object Error]" || p instanceof Error);
    }
    e.isError = re, e.types.isNativeError = re;
    function K(p) {
      return typeof p == "function";
    }
    e.isFunction = K;
    function W(p) {
      return p === null || typeof p == "boolean" || typeof p == "number" || typeof p == "string" || typeof p == "symbol" || // ES6 symbol
      typeof p > "u";
    }
    e.isPrimitive = W, e.isBuffer = kt();
    function $(p) {
      return Object.prototype.toString.call(p);
    }
    function G(p) {
      return p < 10 ? "0" + p.toString(10) : p.toString(10);
    }
    var Y = [
      "Jan",
      "Feb",
      "Mar",
      "Apr",
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Dec"
    ];
    function Z() {
      var p = /* @__PURE__ */ new Date(), m = [
        G(p.getHours()),
        G(p.getMinutes()),
        G(p.getSeconds())
      ].join(":");
      return [p.getDate(), Y[p.getMonth()], m].join(" ");
    }
    e.log = function() {
      console.log("%s - %s", Z(), e.format.apply(e, arguments));
    }, e.inherits = ae, e._extend = function(p, m) {
      if (!m || !C(m))
        return p;
      for (var R = Object.keys(m), F = R.length; F--; )
        p[R[F]] = m[R[F]];
      return p;
    };
    function X(p, m) {
      return Object.prototype.hasOwnProperty.call(p, m);
    }
    var H = typeof Symbol < "u" ? Symbol("util.promisify.custom") : void 0;
    e.promisify = function(m) {
      if (typeof m != "function")
        throw new TypeError('The "original" argument must be of type Function');
      if (H && m[H]) {
        var R = m[H];
        if (typeof R != "function")
          throw new TypeError('The "util.promisify.custom" argument must be of type Function');
        return Object.defineProperty(R, H, {
          value: R,
          enumerable: !1,
          writable: !1,
          configurable: !0
        }), R;
      }
      function R() {
        for (var F, c, u = new Promise(function(B, r) {
          F = B, c = r;
        }), w = [], M = 0; M < arguments.length; M++)
          w.push(arguments[M]);
        w.push(function(B, r) {
          B ? c(B) : F(r);
        });
        try {
          m.apply(this, w);
        } catch (B) {
          c(B);
        }
        return u;
      }
      return Object.setPrototypeOf(R, Object.getPrototypeOf(m)), H && Object.defineProperty(R, H, {
        value: R,
        enumerable: !1,
        writable: !1,
        configurable: !0
      }), Object.defineProperties(
        R,
        n(m)
      );
    }, e.promisify.custom = H;
    function ee(p, m) {
      if (!p) {
        var R = new Error("Promise was rejected with a falsy value");
        R.reason = p, p = R;
      }
      return m(p);
    }
    function te(p) {
      if (typeof p != "function")
        throw new TypeError('The "original" argument must be of type Function');
      function m() {
        for (var R = [], F = 0; F < arguments.length; F++)
          R.push(arguments[F]);
        var c = R.pop();
        if (typeof c != "function")
          throw new TypeError("The last argument must be of type Function");
        var u = this, w = function() {
          return c.apply(u, arguments);
        };
        p.apply(this, R).then(
          function(M) {
            q.nextTick(w.bind(null, null, M));
          },
          function(M) {
            q.nextTick(ee.bind(null, M, w));
          }
        );
      }
      return Object.setPrototypeOf(m, Object.getPrototypeOf(p)), Object.defineProperties(
        m,
        n(p)
      ), m;
    }
    e.callbackify = te;
  }(Ee)), Ee;
}
var He, kr;
function qt() {
  if (kr)
    return He;
  kr = 1;
  function e(_, j) {
    var E = Object.keys(_);
    if (Object.getOwnPropertySymbols) {
      var a = Object.getOwnPropertySymbols(_);
      j && (a = a.filter(function(g) {
        return Object.getOwnPropertyDescriptor(_, g).enumerable;
      })), E.push.apply(E, a);
    }
    return E;
  }
  function n(_) {
    for (var j = 1; j < arguments.length; j++) {
      var E = arguments[j] != null ? arguments[j] : {};
      j % 2 ? e(Object(E), !0).forEach(function(a) {
        i(_, a, E[a]);
      }) : Object.getOwnPropertyDescriptors ? Object.defineProperties(_, Object.getOwnPropertyDescriptors(E)) : e(Object(E)).forEach(function(a) {
        Object.defineProperty(_, a, Object.getOwnPropertyDescriptor(E, a));
      });
    }
    return _;
  }
  function i(_, j, E) {
    return j = f(j), j in _ ? Object.defineProperty(_, j, { value: E, enumerable: !0, configurable: !0, writable: !0 }) : _[j] = E, _;
  }
  function d(_, j) {
    if (!(_ instanceof j))
      throw new TypeError("Cannot call a class as a function");
  }
  function y(_, j) {
    for (var E = 0; E < j.length; E++) {
      var a = j[E];
      a.enumerable = a.enumerable || !1, a.configurable = !0, "value" in a && (a.writable = !0), Object.defineProperty(_, f(a.key), a);
    }
  }
  function s(_, j, E) {
    return j && y(_.prototype, j), E && y(_, E), Object.defineProperty(_, "prototype", { writable: !1 }), _;
  }
  function f(_) {
    var j = l(_, "string");
    return typeof j == "symbol" ? j : String(j);
  }
  function l(_, j) {
    if (typeof _ != "object" || _ === null)
      return _;
    var E = _[Symbol.toPrimitive];
    if (E !== void 0) {
      var a = E.call(_, j || "default");
      if (typeof a != "object")
        return a;
      throw new TypeError("@@toPrimitive must return a primitive value.");
    }
    return (j === "string" ? String : Number)(_);
  }
  var o = ge, h = o.Buffer, v = gt(), S = v.inspect, O = S && S.custom || "inspect";
  function L(_, j, E) {
    h.prototype.copy.call(_, j, E);
  }
  return He = /* @__PURE__ */ function() {
    function _() {
      d(this, _), this.head = null, this.tail = null, this.length = 0;
    }
    return s(_, [{
      key: "push",
      value: function(E) {
        var a = {
          data: E,
          next: null
        };
        this.length > 0 ? this.tail.next = a : this.head = a, this.tail = a, ++this.length;
      }
    }, {
      key: "unshift",
      value: function(E) {
        var a = {
          data: E,
          next: this.head
        };
        this.length === 0 && (this.tail = a), this.head = a, ++this.length;
      }
    }, {
      key: "shift",
      value: function() {
        if (this.length !== 0) {
          var E = this.head.data;
          return this.length === 1 ? this.head = this.tail = null : this.head = this.head.next, --this.length, E;
        }
      }
    }, {
      key: "clear",
      value: function() {
        this.head = this.tail = null, this.length = 0;
      }
    }, {
      key: "join",
      value: function(E) {
        if (this.length === 0)
          return "";
        for (var a = this.head, g = "" + a.data; a = a.next; )
          g += E + a.data;
        return g;
      }
    }, {
      key: "concat",
      value: function(E) {
        if (this.length === 0)
          return h.alloc(0);
        for (var a = h.allocUnsafe(E >>> 0), g = this.head, b = 0; g; )
          L(g.data, a, b), b += g.data.length, g = g.next;
        return a;
      }
      // Consumes a specified amount of bytes or characters from the buffered data.
    }, {
      key: "consume",
      value: function(E, a) {
        var g;
        return E < this.head.data.length ? (g = this.head.data.slice(0, E), this.head.data = this.head.data.slice(E)) : E === this.head.data.length ? g = this.shift() : g = a ? this._getString(E) : this._getBuffer(E), g;
      }
    }, {
      key: "first",
      value: function() {
        return this.head.data;
      }
      // Consumes a specified amount of characters from the buffered data.
    }, {
      key: "_getString",
      value: function(E) {
        var a = this.head, g = 1, b = a.data;
        for (E -= b.length; a = a.next; ) {
          var P = a.data, I = E > P.length ? P.length : E;
          if (I === P.length ? b += P : b += P.slice(0, E), E -= I, E === 0) {
            I === P.length ? (++g, a.next ? this.head = a.next : this.head = this.tail = null) : (this.head = a, a.data = P.slice(I));
            break;
          }
          ++g;
        }
        return this.length -= g, b;
      }
      // Consumes a specified amount of bytes from the buffered data.
    }, {
      key: "_getBuffer",
      value: function(E) {
        var a = h.allocUnsafe(E), g = this.head, b = 1;
        for (g.data.copy(a), E -= g.data.length; g = g.next; ) {
          var P = g.data, I = E > P.length ? P.length : E;
          if (P.copy(a, a.length - E, 0, I), E -= I, E === 0) {
            I === P.length ? (++b, g.next ? this.head = g.next : this.head = this.tail = null) : (this.head = g, g.data = P.slice(I));
            break;
          }
          ++b;
        }
        return this.length -= b, a;
      }
      // Make sure the linked list only shows the minimal necessary information.
    }, {
      key: O,
      value: function(E, a) {
        return S(this, n(n({}, a), {}, {
          // Only inspect one level.
          depth: 0,
          // It should not recurse.
          customInspect: !1
        }));
      }
    }]), _;
  }(), He;
}
var $e, qr;
function vt() {
  if (qr)
    return $e;
  qr = 1;
  function e(f, l) {
    var o = this, h = this._readableState && this._readableState.destroyed, v = this._writableState && this._writableState.destroyed;
    return h || v ? (l ? l(f) : f && (this._writableState ? this._writableState.errorEmitted || (this._writableState.errorEmitted = !0, q.nextTick(y, this, f)) : q.nextTick(y, this, f)), this) : (this._readableState && (this._readableState.destroyed = !0), this._writableState && (this._writableState.destroyed = !0), this._destroy(f || null, function(S) {
      !l && S ? o._writableState ? o._writableState.errorEmitted ? q.nextTick(i, o) : (o._writableState.errorEmitted = !0, q.nextTick(n, o, S)) : q.nextTick(n, o, S) : l ? (q.nextTick(i, o), l(S)) : q.nextTick(i, o);
    }), this);
  }
  function n(f, l) {
    y(f, l), i(f);
  }
  function i(f) {
    f._writableState && !f._writableState.emitClose || f._readableState && !f._readableState.emitClose || f.emit("close");
  }
  function d() {
    this._readableState && (this._readableState.destroyed = !1, this._readableState.reading = !1, this._readableState.ended = !1, this._readableState.endEmitted = !1), this._writableState && (this._writableState.destroyed = !1, this._writableState.ended = !1, this._writableState.ending = !1, this._writableState.finalCalled = !1, this._writableState.prefinished = !1, this._writableState.finished = !1, this._writableState.errorEmitted = !1);
  }
  function y(f, l) {
    f.emit("error", l);
  }
  function s(f, l) {
    var o = f._readableState, h = f._writableState;
    o && o.autoDestroy || h && h.autoDestroy ? f.destroy(l) : f.emit("error", l);
  }
  return $e = {
    destroy: e,
    undestroy: d,
    errorOrDestroy: s
  }, $e;
}
var Ve = {}, Wr;
function ce() {
  if (Wr)
    return Ve;
  Wr = 1;
  function e(l, o) {
    l.prototype = Object.create(o.prototype), l.prototype.constructor = l, l.__proto__ = o;
  }
  var n = {};
  function i(l, o, h) {
    h || (h = Error);
    function v(O, L, _) {
      return typeof o == "string" ? o : o(O, L, _);
    }
    var S = /* @__PURE__ */ function(O) {
      e(L, O);
      function L(_, j, E) {
        return O.call(this, v(_, j, E)) || this;
      }
      return L;
    }(h);
    S.prototype.name = h.name, S.prototype.code = l, n[l] = S;
  }
  function d(l, o) {
    if (Array.isArray(l)) {
      var h = l.length;
      return l = l.map(function(v) {
        return String(v);
      }), h > 2 ? "one of ".concat(o, " ").concat(l.slice(0, h - 1).join(", "), ", or ") + l[h - 1] : h === 2 ? "one of ".concat(o, " ").concat(l[0], " or ").concat(l[1]) : "of ".concat(o, " ").concat(l[0]);
    } else
      return "of ".concat(o, " ").concat(String(l));
  }
  function y(l, o, h) {
    return l.substr(!h || h < 0 ? 0 : +h, o.length) === o;
  }
  function s(l, o, h) {
    return (h === void 0 || h > l.length) && (h = l.length), l.substring(h - o.length, h) === o;
  }
  function f(l, o, h) {
    return typeof h != "number" && (h = 0), h + o.length > l.length ? !1 : l.indexOf(o, h) !== -1;
  }
  return i("ERR_INVALID_OPT_VALUE", function(l, o) {
    return 'The value "' + o + '" is invalid for option "' + l + '"';
  }, TypeError), i("ERR_INVALID_ARG_TYPE", function(l, o, h) {
    var v;
    typeof o == "string" && y(o, "not ") ? (v = "must not be", o = o.replace(/^not /, "")) : v = "must be";
    var S;
    if (s(l, " argument"))
      S = "The ".concat(l, " ").concat(v, " ").concat(d(o, "type"));
    else {
      var O = f(l, ".") ? "property" : "argument";
      S = 'The "'.concat(l, '" ').concat(O, " ").concat(v, " ").concat(d(o, "type"));
    }
    return S += ". Received type ".concat(typeof h), S;
  }, TypeError), i("ERR_STREAM_PUSH_AFTER_EOF", "stream.push() after EOF"), i("ERR_METHOD_NOT_IMPLEMENTED", function(l) {
    return "The " + l + " method is not implemented";
  }), i("ERR_STREAM_PREMATURE_CLOSE", "Premature close"), i("ERR_STREAM_DESTROYED", function(l) {
    return "Cannot call " + l + " after a stream was destroyed";
  }), i("ERR_MULTIPLE_CALLBACK", "Callback called multiple times"), i("ERR_STREAM_CANNOT_PIPE", "Cannot pipe, not readable"), i("ERR_STREAM_WRITE_AFTER_END", "write after end"), i("ERR_STREAM_NULL_VALUES", "May not write null values to stream", TypeError), i("ERR_UNKNOWN_ENCODING", function(l) {
    return "Unknown encoding: " + l;
  }, TypeError), i("ERR_STREAM_UNSHIFT_AFTER_END_EVENT", "stream.unshift() after end event"), Ve.codes = n, Ve;
}
var ze, xr;
function bt() {
  if (xr)
    return ze;
  xr = 1;
  var e = ce().codes.ERR_INVALID_OPT_VALUE;
  function n(d, y, s) {
    return d.highWaterMark != null ? d.highWaterMark : y ? d[s] : null;
  }
  function i(d, y, s, f) {
    var l = n(y, f, s);
    if (l != null) {
      if (!(isFinite(l) && Math.floor(l) === l) || l < 0) {
        var o = f ? s : "highWaterMark";
        throw new e(o, l);
      }
      return Math.floor(l);
    }
    return d.objectMode ? 16 : 16 * 1024;
  }
  return ze = {
    getHighWaterMark: i
  }, ze;
}
var Ke, Gr;
function Wt() {
  if (Gr)
    return Ke;
  Gr = 1, Ke = e;
  function e(i, d) {
    if (n("noDeprecation"))
      return i;
    var y = !1;
    function s() {
      if (!y) {
        if (n("throwDeprecation"))
          throw new Error(d);
        n("traceDeprecation") ? console.trace(d) : console.warn(d), y = !0;
      }
      return i.apply(this, arguments);
    }
    return s;
  }
  function n(i) {
    try {
      if (!ie.localStorage)
        return !1;
    } catch {
      return !1;
    }
    var d = ie.localStorage[i];
    return d == null ? !1 : String(d).toLowerCase() === "true";
  }
  return Ke;
}
var Je, Hr;
function mt() {
  if (Hr)
    return Je;
  Hr = 1, Je = U;
  function e(c) {
    var u = this;
    this.next = null, this.entry = null, this.finish = function() {
      F(u, c);
    };
  }
  var n;
  U.WritableState = D;
  var i = {
    deprecate: Wt()
  }, d = lt(), y = ge.Buffer, s = (typeof ie < "u" ? ie : typeof window < "u" ? window : typeof self < "u" ? self : {}).Uint8Array || function() {
  };
  function f(c) {
    return y.from(c);
  }
  function l(c) {
    return y.isBuffer(c) || c instanceof s;
  }
  var o = vt(), h = bt(), v = h.getHighWaterMark, S = ce().codes, O = S.ERR_INVALID_ARG_TYPE, L = S.ERR_METHOD_NOT_IMPLEMENTED, _ = S.ERR_MULTIPLE_CALLBACK, j = S.ERR_STREAM_CANNOT_PIPE, E = S.ERR_STREAM_DESTROYED, a = S.ERR_STREAM_NULL_VALUES, g = S.ERR_STREAM_WRITE_AFTER_END, b = S.ERR_UNKNOWN_ENCODING, P = o.errorOrDestroy;
  ae(U, d);
  function I() {
  }
  function D(c, u, w) {
    n = n || le(), c = c || {}, typeof w != "boolean" && (w = u instanceof n), this.objectMode = !!c.objectMode, w && (this.objectMode = this.objectMode || !!c.writableObjectMode), this.highWaterMark = v(this, c, "writableHighWaterMark", w), this.finalCalled = !1, this.needDrain = !1, this.ending = !1, this.ended = !1, this.finished = !1, this.destroyed = !1;
    var M = c.decodeStrings === !1;
    this.decodeStrings = !M, this.defaultEncoding = c.defaultEncoding || "utf8", this.length = 0, this.writing = !1, this.corked = 0, this.sync = !0, this.bufferProcessing = !1, this.onwrite = function(B) {
      Y(u, B);
    }, this.writecb = null, this.writelen = 0, this.bufferedRequest = null, this.lastBufferedRequest = null, this.pendingcb = 0, this.prefinished = !1, this.errorEmitted = !1, this.emitClose = c.emitClose !== !1, this.autoDestroy = !!c.autoDestroy, this.bufferedRequestCount = 0, this.corkedRequestsFree = new e(this);
  }
  D.prototype.getBuffer = function() {
    for (var u = this.bufferedRequest, w = []; u; )
      w.push(u), u = u.next;
    return w;
  }, function() {
    try {
      Object.defineProperty(D.prototype, "buffer", {
        get: i.deprecate(function() {
          return this.getBuffer();
        }, "_writableState.buffer is deprecated. Use _writableState.getBuffer instead.", "DEP0003")
      });
    } catch {
    }
  }();
  var x;
  typeof Symbol == "function" && Symbol.hasInstance && typeof Function.prototype[Symbol.hasInstance] == "function" ? (x = Function.prototype[Symbol.hasInstance], Object.defineProperty(U, Symbol.hasInstance, {
    value: function(u) {
      return x.call(this, u) ? !0 : this !== U ? !1 : u && u._writableState instanceof D;
    }
  })) : x = function(u) {
    return u instanceof this;
  };
  function U(c) {
    n = n || le();
    var u = this instanceof n;
    if (!u && !x.call(U, this))
      return new U(c);
    this._writableState = new D(c, this, u), this.writable = !0, c && (typeof c.write == "function" && (this._write = c.write), typeof c.writev == "function" && (this._writev = c.writev), typeof c.destroy == "function" && (this._destroy = c.destroy), typeof c.final == "function" && (this._final = c.final)), d.call(this);
  }
  U.prototype.pipe = function() {
    P(this, new j());
  };
  function C(c, u) {
    var w = new g();
    P(c, w), q.nextTick(u, w);
  }
  function Q(c, u, w, M) {
    var B;
    return w === null ? B = new a() : typeof w != "string" && !u.objectMode && (B = new O("chunk", ["string", "Buffer"], w)), B ? (P(c, B), q.nextTick(M, B), !1) : !0;
  }
  U.prototype.write = function(c, u, w) {
    var M = this._writableState, B = !1, r = !M.objectMode && l(c);
    return r && !y.isBuffer(c) && (c = f(c)), typeof u == "function" && (w = u, u = null), r ? u = "buffer" : u || (u = M.defaultEncoding), typeof w != "function" && (w = I), M.ending ? C(this, w) : (r || Q(this, M, c, w)) && (M.pendingcb++, B = K(this, M, r, c, u, w)), B;
  }, U.prototype.cork = function() {
    this._writableState.corked++;
  }, U.prototype.uncork = function() {
    var c = this._writableState;
    c.corked && (c.corked--, !c.writing && !c.corked && !c.bufferProcessing && c.bufferedRequest && H(this, c));
  }, U.prototype.setDefaultEncoding = function(u) {
    if (typeof u == "string" && (u = u.toLowerCase()), !(["hex", "utf8", "utf-8", "ascii", "binary", "base64", "ucs2", "ucs-2", "utf16le", "utf-16le", "raw"].indexOf((u + "").toLowerCase()) > -1))
      throw new b(u);
    return this._writableState.defaultEncoding = u, this;
  }, Object.defineProperty(U.prototype, "writableBuffer", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState && this._writableState.getBuffer();
    }
  });
  function re(c, u, w) {
    return !c.objectMode && c.decodeStrings !== !1 && typeof u == "string" && (u = y.from(u, w)), u;
  }
  Object.defineProperty(U.prototype, "writableHighWaterMark", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState.highWaterMark;
    }
  });
  function K(c, u, w, M, B, r) {
    if (!w) {
      var t = re(u, M, B);
      M !== t && (w = !0, B = "buffer", M = t);
    }
    var T = u.objectMode ? 1 : M.length;
    u.length += T;
    var N = u.length < u.highWaterMark;
    if (N || (u.needDrain = !0), u.writing || u.corked) {
      var J = u.lastBufferedRequest;
      u.lastBufferedRequest = {
        chunk: M,
        encoding: B,
        isBuf: w,
        callback: r,
        next: null
      }, J ? J.next = u.lastBufferedRequest : u.bufferedRequest = u.lastBufferedRequest, u.bufferedRequestCount += 1;
    } else
      W(c, u, !1, T, M, B, r);
    return N;
  }
  function W(c, u, w, M, B, r, t) {
    u.writelen = M, u.writecb = t, u.writing = !0, u.sync = !0, u.destroyed ? u.onwrite(new E("write")) : w ? c._writev(B, u.onwrite) : c._write(B, r, u.onwrite), u.sync = !1;
  }
  function $(c, u, w, M, B) {
    --u.pendingcb, w ? (q.nextTick(B, M), q.nextTick(m, c, u), c._writableState.errorEmitted = !0, P(c, M)) : (B(M), c._writableState.errorEmitted = !0, P(c, M), m(c, u));
  }
  function G(c) {
    c.writing = !1, c.writecb = null, c.length -= c.writelen, c.writelen = 0;
  }
  function Y(c, u) {
    var w = c._writableState, M = w.sync, B = w.writecb;
    if (typeof B != "function")
      throw new _();
    if (G(w), u)
      $(c, w, M, u, B);
    else {
      var r = ee(w) || c.destroyed;
      !r && !w.corked && !w.bufferProcessing && w.bufferedRequest && H(c, w), M ? q.nextTick(Z, c, w, r, B) : Z(c, w, r, B);
    }
  }
  function Z(c, u, w, M) {
    w || X(c, u), u.pendingcb--, M(), m(c, u);
  }
  function X(c, u) {
    u.length === 0 && u.needDrain && (u.needDrain = !1, c.emit("drain"));
  }
  function H(c, u) {
    u.bufferProcessing = !0;
    var w = u.bufferedRequest;
    if (c._writev && w && w.next) {
      var M = u.bufferedRequestCount, B = new Array(M), r = u.corkedRequestsFree;
      r.entry = w;
      for (var t = 0, T = !0; w; )
        B[t] = w, w.isBuf || (T = !1), w = w.next, t += 1;
      B.allBuffers = T, W(c, u, !0, u.length, B, "", r.finish), u.pendingcb++, u.lastBufferedRequest = null, r.next ? (u.corkedRequestsFree = r.next, r.next = null) : u.corkedRequestsFree = new e(u), u.bufferedRequestCount = 0;
    } else {
      for (; w; ) {
        var N = w.chunk, J = w.encoding, k = w.callback, V = u.objectMode ? 1 : N.length;
        if (W(c, u, !1, V, N, J, k), w = w.next, u.bufferedRequestCount--, u.writing)
          break;
      }
      w === null && (u.lastBufferedRequest = null);
    }
    u.bufferedRequest = w, u.bufferProcessing = !1;
  }
  U.prototype._write = function(c, u, w) {
    w(new L("_write()"));
  }, U.prototype._writev = null, U.prototype.end = function(c, u, w) {
    var M = this._writableState;
    return typeof c == "function" ? (w = c, c = null, u = null) : typeof u == "function" && (w = u, u = null), c != null && this.write(c, u), M.corked && (M.corked = 1, this.uncork()), M.ending || R(this, M, w), this;
  }, Object.defineProperty(U.prototype, "writableLength", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState.length;
    }
  });
  function ee(c) {
    return c.ending && c.length === 0 && c.bufferedRequest === null && !c.finished && !c.writing;
  }
  function te(c, u) {
    c._final(function(w) {
      u.pendingcb--, w && P(c, w), u.prefinished = !0, c.emit("prefinish"), m(c, u);
    });
  }
  function p(c, u) {
    !u.prefinished && !u.finalCalled && (typeof c._final == "function" && !u.destroyed ? (u.pendingcb++, u.finalCalled = !0, q.nextTick(te, c, u)) : (u.prefinished = !0, c.emit("prefinish")));
  }
  function m(c, u) {
    var w = ee(u);
    if (w && (p(c, u), u.pendingcb === 0 && (u.finished = !0, c.emit("finish"), u.autoDestroy))) {
      var M = c._readableState;
      (!M || M.autoDestroy && M.endEmitted) && c.destroy();
    }
    return w;
  }
  function R(c, u, w) {
    u.ending = !0, m(c, u), w && (u.finished ? q.nextTick(w) : c.once("finish", w)), u.ended = !0, c.writable = !1;
  }
  function F(c, u, w) {
    var M = c.entry;
    for (c.entry = null; M; ) {
      var B = M.callback;
      u.pendingcb--, B(w), M = M.next;
    }
    u.corkedRequestsFree.next = c;
  }
  return Object.defineProperty(U.prototype, "destroyed", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState === void 0 ? !1 : this._writableState.destroyed;
    },
    set: function(u) {
      this._writableState && (this._writableState.destroyed = u);
    }
  }), U.prototype.destroy = o.destroy, U.prototype._undestroy = o.undestroy, U.prototype._destroy = function(c, u) {
    u(c);
  }, Je;
}
var Ye, $r;
function le() {
  if ($r)
    return Ye;
  $r = 1;
  var e = Object.keys || function(h) {
    var v = [];
    for (var S in h)
      v.push(S);
    return v;
  };
  Ye = f;
  var n = wt(), i = mt();
  ae(f, n);
  for (var d = e(i.prototype), y = 0; y < d.length; y++) {
    var s = d[y];
    f.prototype[s] || (f.prototype[s] = i.prototype[s]);
  }
  function f(h) {
    if (!(this instanceof f))
      return new f(h);
    n.call(this, h), i.call(this, h), this.allowHalfOpen = !0, h && (h.readable === !1 && (this.readable = !1), h.writable === !1 && (this.writable = !1), h.allowHalfOpen === !1 && (this.allowHalfOpen = !1, this.once("end", l)));
  }
  Object.defineProperty(f.prototype, "writableHighWaterMark", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState.highWaterMark;
    }
  }), Object.defineProperty(f.prototype, "writableBuffer", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState && this._writableState.getBuffer();
    }
  }), Object.defineProperty(f.prototype, "writableLength", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._writableState.length;
    }
  });
  function l() {
    this._writableState.ended || q.nextTick(o, this);
  }
  function o(h) {
    h.end();
  }
  return Object.defineProperty(f.prototype, "destroyed", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._readableState === void 0 || this._writableState === void 0 ? !1 : this._readableState.destroyed && this._writableState.destroyed;
    },
    set: function(v) {
      this._readableState === void 0 || this._writableState === void 0 || (this._readableState.destroyed = v, this._writableState.destroyed = v);
    }
  }), Ye;
}
var Ze = {}, ye = { exports: {} };
/*! safe-buffer. MIT License. Feross Aboukhadijeh <https://feross.org/opensource> */
var Vr;
function xt() {
  return Vr || (Vr = 1, function(e, n) {
    var i = ge, d = i.Buffer;
    function y(f, l) {
      for (var o in f)
        l[o] = f[o];
    }
    d.from && d.alloc && d.allocUnsafe && d.allocUnsafeSlow ? e.exports = i : (y(i, n), n.Buffer = s);
    function s(f, l, o) {
      return d(f, l, o);
    }
    s.prototype = Object.create(d.prototype), y(d, s), s.from = function(f, l, o) {
      if (typeof f == "number")
        throw new TypeError("Argument must not be a number");
      return d(f, l, o);
    }, s.alloc = function(f, l, o) {
      if (typeof f != "number")
        throw new TypeError("Argument must be a number");
      var h = d(f);
      return l !== void 0 ? typeof o == "string" ? h.fill(l, o) : h.fill(l) : h.fill(0), h;
    }, s.allocUnsafe = function(f) {
      if (typeof f != "number")
        throw new TypeError("Argument must be a number");
      return d(f);
    }, s.allocUnsafeSlow = function(f) {
      if (typeof f != "number")
        throw new TypeError("Argument must be a number");
      return i.SlowBuffer(f);
    };
  }(ye, ye.exports)), ye.exports;
}
var zr;
function Kr() {
  if (zr)
    return Ze;
  zr = 1;
  var e = xt().Buffer, n = e.isEncoding || function(a) {
    switch (a = "" + a, a && a.toLowerCase()) {
      case "hex":
      case "utf8":
      case "utf-8":
      case "ascii":
      case "binary":
      case "base64":
      case "ucs2":
      case "ucs-2":
      case "utf16le":
      case "utf-16le":
      case "raw":
        return !0;
      default:
        return !1;
    }
  };
  function i(a) {
    if (!a)
      return "utf8";
    for (var g; ; )
      switch (a) {
        case "utf8":
        case "utf-8":
          return "utf8";
        case "ucs2":
        case "ucs-2":
        case "utf16le":
        case "utf-16le":
          return "utf16le";
        case "latin1":
        case "binary":
          return "latin1";
        case "base64":
        case "ascii":
        case "hex":
          return a;
        default:
          if (g)
            return;
          a = ("" + a).toLowerCase(), g = !0;
      }
  }
  function d(a) {
    var g = i(a);
    if (typeof g != "string" && (e.isEncoding === n || !n(a)))
      throw new Error("Unknown encoding: " + a);
    return g || a;
  }
  Ze.StringDecoder = y;
  function y(a) {
    this.encoding = d(a);
    var g;
    switch (this.encoding) {
      case "utf16le":
        this.text = S, this.end = O, g = 4;
        break;
      case "utf8":
        this.fillLast = o, g = 4;
        break;
      case "base64":
        this.text = L, this.end = _, g = 3;
        break;
      default:
        this.write = j, this.end = E;
        return;
    }
    this.lastNeed = 0, this.lastTotal = 0, this.lastChar = e.allocUnsafe(g);
  }
  y.prototype.write = function(a) {
    if (a.length === 0)
      return "";
    var g, b;
    if (this.lastNeed) {
      if (g = this.fillLast(a), g === void 0)
        return "";
      b = this.lastNeed, this.lastNeed = 0;
    } else
      b = 0;
    return b < a.length ? g ? g + this.text(a, b) : this.text(a, b) : g || "";
  }, y.prototype.end = v, y.prototype.text = h, y.prototype.fillLast = function(a) {
    if (this.lastNeed <= a.length)
      return a.copy(this.lastChar, this.lastTotal - this.lastNeed, 0, this.lastNeed), this.lastChar.toString(this.encoding, 0, this.lastTotal);
    a.copy(this.lastChar, this.lastTotal - this.lastNeed, 0, a.length), this.lastNeed -= a.length;
  };
  function s(a) {
    return a <= 127 ? 0 : a >> 5 === 6 ? 2 : a >> 4 === 14 ? 3 : a >> 3 === 30 ? 4 : a >> 6 === 2 ? -1 : -2;
  }
  function f(a, g, b) {
    var P = g.length - 1;
    if (P < b)
      return 0;
    var I = s(g[P]);
    return I >= 0 ? (I > 0 && (a.lastNeed = I - 1), I) : --P < b || I === -2 ? 0 : (I = s(g[P]), I >= 0 ? (I > 0 && (a.lastNeed = I - 2), I) : --P < b || I === -2 ? 0 : (I = s(g[P]), I >= 0 ? (I > 0 && (I === 2 ? I = 0 : a.lastNeed = I - 3), I) : 0));
  }
  function l(a, g, b) {
    if ((g[0] & 192) !== 128)
      return a.lastNeed = 0, "�";
    if (a.lastNeed > 1 && g.length > 1) {
      if ((g[1] & 192) !== 128)
        return a.lastNeed = 1, "�";
      if (a.lastNeed > 2 && g.length > 2 && (g[2] & 192) !== 128)
        return a.lastNeed = 2, "�";
    }
  }
  function o(a) {
    var g = this.lastTotal - this.lastNeed, b = l(this, a);
    if (b !== void 0)
      return b;
    if (this.lastNeed <= a.length)
      return a.copy(this.lastChar, g, 0, this.lastNeed), this.lastChar.toString(this.encoding, 0, this.lastTotal);
    a.copy(this.lastChar, g, 0, a.length), this.lastNeed -= a.length;
  }
  function h(a, g) {
    var b = f(this, a, g);
    if (!this.lastNeed)
      return a.toString("utf8", g);
    this.lastTotal = b;
    var P = a.length - (b - this.lastNeed);
    return a.copy(this.lastChar, 0, P), a.toString("utf8", g, P);
  }
  function v(a) {
    var g = a && a.length ? this.write(a) : "";
    return this.lastNeed ? g + "�" : g;
  }
  function S(a, g) {
    if ((a.length - g) % 2 === 0) {
      var b = a.toString("utf16le", g);
      if (b) {
        var P = b.charCodeAt(b.length - 1);
        if (P >= 55296 && P <= 56319)
          return this.lastNeed = 2, this.lastTotal = 4, this.lastChar[0] = a[a.length - 2], this.lastChar[1] = a[a.length - 1], b.slice(0, -1);
      }
      return b;
    }
    return this.lastNeed = 1, this.lastTotal = 2, this.lastChar[0] = a[a.length - 1], a.toString("utf16le", g, a.length - 1);
  }
  function O(a) {
    var g = a && a.length ? this.write(a) : "";
    if (this.lastNeed) {
      var b = this.lastTotal - this.lastNeed;
      return g + this.lastChar.toString("utf16le", 0, b);
    }
    return g;
  }
  function L(a, g) {
    var b = (a.length - g) % 3;
    return b === 0 ? a.toString("base64", g) : (this.lastNeed = 3 - b, this.lastTotal = 3, b === 1 ? this.lastChar[0] = a[a.length - 1] : (this.lastChar[0] = a[a.length - 2], this.lastChar[1] = a[a.length - 1]), a.toString("base64", g, a.length - b));
  }
  function _(a) {
    var g = a && a.length ? this.write(a) : "";
    return this.lastNeed ? g + this.lastChar.toString("base64", 0, 3 - this.lastNeed) : g;
  }
  function j(a) {
    return a.toString(this.encoding);
  }
  function E(a) {
    return a && a.length ? this.write(a) : "";
  }
  return Ze;
}
var Qe, Jr;
function cr() {
  if (Jr)
    return Qe;
  Jr = 1;
  var e = ce().codes.ERR_STREAM_PREMATURE_CLOSE;
  function n(s) {
    var f = !1;
    return function() {
      if (!f) {
        f = !0;
        for (var l = arguments.length, o = new Array(l), h = 0; h < l; h++)
          o[h] = arguments[h];
        s.apply(this, o);
      }
    };
  }
  function i() {
  }
  function d(s) {
    return s.setHeader && typeof s.abort == "function";
  }
  function y(s, f, l) {
    if (typeof f == "function")
      return y(s, null, f);
    f || (f = {}), l = n(l || i);
    var o = f.readable || f.readable !== !1 && s.readable, h = f.writable || f.writable !== !1 && s.writable, v = function() {
      s.writable || O();
    }, S = s._writableState && s._writableState.finished, O = function() {
      h = !1, S = !0, o || l.call(s);
    }, L = s._readableState && s._readableState.endEmitted, _ = function() {
      o = !1, L = !0, h || l.call(s);
    }, j = function(b) {
      l.call(s, b);
    }, E = function() {
      var b;
      if (o && !L)
        return (!s._readableState || !s._readableState.ended) && (b = new e()), l.call(s, b);
      if (h && !S)
        return (!s._writableState || !s._writableState.ended) && (b = new e()), l.call(s, b);
    }, a = function() {
      s.req.on("finish", O);
    };
    return d(s) ? (s.on("complete", O), s.on("abort", E), s.req ? a() : s.on("request", a)) : h && !s._writableState && (s.on("end", v), s.on("close", v)), s.on("end", _), s.on("finish", O), f.error !== !1 && s.on("error", j), s.on("close", E), function() {
      s.removeListener("complete", O), s.removeListener("abort", E), s.removeListener("request", a), s.req && s.req.removeListener("finish", O), s.removeListener("end", v), s.removeListener("close", v), s.removeListener("finish", O), s.removeListener("end", _), s.removeListener("error", j), s.removeListener("close", E);
    };
  }
  return Qe = y, Qe;
}
var Xe, Yr;
function Gt() {
  if (Yr)
    return Xe;
  Yr = 1;
  var e;
  function n(b, P, I) {
    return P = i(P), P in b ? Object.defineProperty(b, P, { value: I, enumerable: !0, configurable: !0, writable: !0 }) : b[P] = I, b;
  }
  function i(b) {
    var P = d(b, "string");
    return typeof P == "symbol" ? P : String(P);
  }
  function d(b, P) {
    if (typeof b != "object" || b === null)
      return b;
    var I = b[Symbol.toPrimitive];
    if (I !== void 0) {
      var D = I.call(b, P || "default");
      if (typeof D != "object")
        return D;
      throw new TypeError("@@toPrimitive must return a primitive value.");
    }
    return (P === "string" ? String : Number)(b);
  }
  var y = cr(), s = Symbol("lastResolve"), f = Symbol("lastReject"), l = Symbol("error"), o = Symbol("ended"), h = Symbol("lastPromise"), v = Symbol("handlePromise"), S = Symbol("stream");
  function O(b, P) {
    return {
      value: b,
      done: P
    };
  }
  function L(b) {
    var P = b[s];
    if (P !== null) {
      var I = b[S].read();
      I !== null && (b[h] = null, b[s] = null, b[f] = null, P(O(I, !1)));
    }
  }
  function _(b) {
    q.nextTick(L, b);
  }
  function j(b, P) {
    return function(I, D) {
      b.then(function() {
        if (P[o]) {
          I(O(void 0, !0));
          return;
        }
        P[v](I, D);
      }, D);
    };
  }
  var E = Object.getPrototypeOf(function() {
  }), a = Object.setPrototypeOf((e = {
    get stream() {
      return this[S];
    },
    next: function() {
      var P = this, I = this[l];
      if (I !== null)
        return Promise.reject(I);
      if (this[o])
        return Promise.resolve(O(void 0, !0));
      if (this[S].destroyed)
        return new Promise(function(C, Q) {
          q.nextTick(function() {
            P[l] ? Q(P[l]) : C(O(void 0, !0));
          });
        });
      var D = this[h], x;
      if (D)
        x = new Promise(j(D, this));
      else {
        var U = this[S].read();
        if (U !== null)
          return Promise.resolve(O(U, !1));
        x = new Promise(this[v]);
      }
      return this[h] = x, x;
    }
  }, n(e, Symbol.asyncIterator, function() {
    return this;
  }), n(e, "return", function() {
    var P = this;
    return new Promise(function(I, D) {
      P[S].destroy(null, function(x) {
        if (x) {
          D(x);
          return;
        }
        I(O(void 0, !0));
      });
    });
  }), e), E), g = function(P) {
    var I, D = Object.create(a, (I = {}, n(I, S, {
      value: P,
      writable: !0
    }), n(I, s, {
      value: null,
      writable: !0
    }), n(I, f, {
      value: null,
      writable: !0
    }), n(I, l, {
      value: null,
      writable: !0
    }), n(I, o, {
      value: P._readableState.endEmitted,
      writable: !0
    }), n(I, v, {
      value: function(U, C) {
        var Q = D[S].read();
        Q ? (D[h] = null, D[s] = null, D[f] = null, U(O(Q, !1))) : (D[s] = U, D[f] = C);
      },
      writable: !0
    }), I));
    return D[h] = null, y(P, function(x) {
      if (x && x.code !== "ERR_STREAM_PREMATURE_CLOSE") {
        var U = D[f];
        U !== null && (D[h] = null, D[s] = null, D[f] = null, U(x)), D[l] = x;
        return;
      }
      var C = D[s];
      C !== null && (D[h] = null, D[s] = null, D[f] = null, C(O(void 0, !0))), D[o] = !0;
    }), P.on("readable", _.bind(null, D)), D;
  };
  return Xe = g, Xe;
}
var er, Zr;
function Ht() {
  return Zr || (Zr = 1, er = function() {
    throw new Error("Readable.from is not available in the browser");
  }), er;
}
var rr, Qr;
function wt() {
  if (Qr)
    return rr;
  Qr = 1, rr = C;
  var e;
  C.ReadableState = U, fr.EventEmitter;
  var n = function(t, T) {
    return t.listeners(T).length;
  }, i = lt(), d = ge.Buffer, y = (typeof ie < "u" ? ie : typeof window < "u" ? window : typeof self < "u" ? self : {}).Uint8Array || function() {
  };
  function s(r) {
    return d.from(r);
  }
  function f(r) {
    return d.isBuffer(r) || r instanceof y;
  }
  var l = gt(), o;
  l && l.debuglog ? o = l.debuglog("stream") : o = function() {
  };
  var h = qt(), v = vt(), S = bt(), O = S.getHighWaterMark, L = ce().codes, _ = L.ERR_INVALID_ARG_TYPE, j = L.ERR_STREAM_PUSH_AFTER_EOF, E = L.ERR_METHOD_NOT_IMPLEMENTED, a = L.ERR_STREAM_UNSHIFT_AFTER_END_EVENT, g, b, P;
  ae(C, i);
  var I = v.errorOrDestroy, D = ["error", "close", "destroy", "pause", "resume"];
  function x(r, t, T) {
    if (typeof r.prependListener == "function")
      return r.prependListener(t, T);
    !r._events || !r._events[t] ? r.on(t, T) : Array.isArray(r._events[t]) ? r._events[t].unshift(T) : r._events[t] = [T, r._events[t]];
  }
  function U(r, t, T) {
    e = e || le(), r = r || {}, typeof T != "boolean" && (T = t instanceof e), this.objectMode = !!r.objectMode, T && (this.objectMode = this.objectMode || !!r.readableObjectMode), this.highWaterMark = O(this, r, "readableHighWaterMark", T), this.buffer = new h(), this.length = 0, this.pipes = null, this.pipesCount = 0, this.flowing = null, this.ended = !1, this.endEmitted = !1, this.reading = !1, this.sync = !0, this.needReadable = !1, this.emittedReadable = !1, this.readableListening = !1, this.resumeScheduled = !1, this.paused = !0, this.emitClose = r.emitClose !== !1, this.autoDestroy = !!r.autoDestroy, this.destroyed = !1, this.defaultEncoding = r.defaultEncoding || "utf8", this.awaitDrain = 0, this.readingMore = !1, this.decoder = null, this.encoding = null, r.encoding && (g || (g = Kr().StringDecoder), this.decoder = new g(r.encoding), this.encoding = r.encoding);
  }
  function C(r) {
    if (e = e || le(), !(this instanceof C))
      return new C(r);
    var t = this instanceof e;
    this._readableState = new U(r, this, t), this.readable = !0, r && (typeof r.read == "function" && (this._read = r.read), typeof r.destroy == "function" && (this._destroy = r.destroy)), i.call(this);
  }
  Object.defineProperty(C.prototype, "destroyed", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._readableState === void 0 ? !1 : this._readableState.destroyed;
    },
    set: function(t) {
      this._readableState && (this._readableState.destroyed = t);
    }
  }), C.prototype.destroy = v.destroy, C.prototype._undestroy = v.undestroy, C.prototype._destroy = function(r, t) {
    t(r);
  }, C.prototype.push = function(r, t) {
    var T = this._readableState, N;
    return T.objectMode ? N = !0 : typeof r == "string" && (t = t || T.defaultEncoding, t !== T.encoding && (r = d.from(r, t), t = ""), N = !0), Q(this, r, t, !1, N);
  }, C.prototype.unshift = function(r) {
    return Q(this, r, null, !0, !1);
  };
  function Q(r, t, T, N, J) {
    o("readableAddChunk", t);
    var k = r._readableState;
    if (t === null)
      k.reading = !1, Y(r, k);
    else {
      var V;
      if (J || (V = K(k, t)), V)
        I(r, V);
      else if (k.objectMode || t && t.length > 0)
        if (typeof t != "string" && !k.objectMode && Object.getPrototypeOf(t) !== d.prototype && (t = s(t)), N)
          k.endEmitted ? I(r, new a()) : re(r, k, t, !0);
        else if (k.ended)
          I(r, new j());
        else {
          if (k.destroyed)
            return !1;
          k.reading = !1, k.decoder && !T ? (t = k.decoder.write(t), k.objectMode || t.length !== 0 ? re(r, k, t, !1) : H(r, k)) : re(r, k, t, !1);
        }
      else
        N || (k.reading = !1, H(r, k));
    }
    return !k.ended && (k.length < k.highWaterMark || k.length === 0);
  }
  function re(r, t, T, N) {
    t.flowing && t.length === 0 && !t.sync ? (t.awaitDrain = 0, r.emit("data", T)) : (t.length += t.objectMode ? 1 : T.length, N ? t.buffer.unshift(T) : t.buffer.push(T), t.needReadable && Z(r)), H(r, t);
  }
  function K(r, t) {
    var T;
    return !f(t) && typeof t != "string" && t !== void 0 && !r.objectMode && (T = new _("chunk", ["string", "Buffer", "Uint8Array"], t)), T;
  }
  C.prototype.isPaused = function() {
    return this._readableState.flowing === !1;
  }, C.prototype.setEncoding = function(r) {
    g || (g = Kr().StringDecoder);
    var t = new g(r);
    this._readableState.decoder = t, this._readableState.encoding = this._readableState.decoder.encoding;
    for (var T = this._readableState.buffer.head, N = ""; T !== null; )
      N += t.write(T.data), T = T.next;
    return this._readableState.buffer.clear(), N !== "" && this._readableState.buffer.push(N), this._readableState.length = N.length, this;
  };
  var W = 1073741824;
  function $(r) {
    return r >= W ? r = W : (r--, r |= r >>> 1, r |= r >>> 2, r |= r >>> 4, r |= r >>> 8, r |= r >>> 16, r++), r;
  }
  function G(r, t) {
    return r <= 0 || t.length === 0 && t.ended ? 0 : t.objectMode ? 1 : r !== r ? t.flowing && t.length ? t.buffer.head.data.length : t.length : (r > t.highWaterMark && (t.highWaterMark = $(r)), r <= t.length ? r : t.ended ? t.length : (t.needReadable = !0, 0));
  }
  C.prototype.read = function(r) {
    o("read", r), r = parseInt(r, 10);
    var t = this._readableState, T = r;
    if (r !== 0 && (t.emittedReadable = !1), r === 0 && t.needReadable && ((t.highWaterMark !== 0 ? t.length >= t.highWaterMark : t.length > 0) || t.ended))
      return o("read: emitReadable", t.length, t.ended), t.length === 0 && t.ended ? w(this) : Z(this), null;
    if (r = G(r, t), r === 0 && t.ended)
      return t.length === 0 && w(this), null;
    var N = t.needReadable;
    o("need readable", N), (t.length === 0 || t.length - r < t.highWaterMark) && (N = !0, o("length less than watermark", N)), t.ended || t.reading ? (N = !1, o("reading or ended", N)) : N && (o("do read"), t.reading = !0, t.sync = !0, t.length === 0 && (t.needReadable = !0), this._read(t.highWaterMark), t.sync = !1, t.reading || (r = G(T, t)));
    var J;
    return r > 0 ? J = u(r, t) : J = null, J === null ? (t.needReadable = t.length <= t.highWaterMark, r = 0) : (t.length -= r, t.awaitDrain = 0), t.length === 0 && (t.ended || (t.needReadable = !0), T !== r && t.ended && w(this)), J !== null && this.emit("data", J), J;
  };
  function Y(r, t) {
    if (o("onEofChunk"), !t.ended) {
      if (t.decoder) {
        var T = t.decoder.end();
        T && T.length && (t.buffer.push(T), t.length += t.objectMode ? 1 : T.length);
      }
      t.ended = !0, t.sync ? Z(r) : (t.needReadable = !1, t.emittedReadable || (t.emittedReadable = !0, X(r)));
    }
  }
  function Z(r) {
    var t = r._readableState;
    o("emitReadable", t.needReadable, t.emittedReadable), t.needReadable = !1, t.emittedReadable || (o("emitReadable", t.flowing), t.emittedReadable = !0, q.nextTick(X, r));
  }
  function X(r) {
    var t = r._readableState;
    o("emitReadable_", t.destroyed, t.length, t.ended), !t.destroyed && (t.length || t.ended) && (r.emit("readable"), t.emittedReadable = !1), t.needReadable = !t.flowing && !t.ended && t.length <= t.highWaterMark, c(r);
  }
  function H(r, t) {
    t.readingMore || (t.readingMore = !0, q.nextTick(ee, r, t));
  }
  function ee(r, t) {
    for (; !t.reading && !t.ended && (t.length < t.highWaterMark || t.flowing && t.length === 0); ) {
      var T = t.length;
      if (o("maybeReadMore read 0"), r.read(0), T === t.length)
        break;
    }
    t.readingMore = !1;
  }
  C.prototype._read = function(r) {
    I(this, new E("_read()"));
  }, C.prototype.pipe = function(r, t) {
    var T = this, N = this._readableState;
    switch (N.pipesCount) {
      case 0:
        N.pipes = r;
        break;
      case 1:
        N.pipes = [N.pipes, r];
        break;
      default:
        N.pipes.push(r);
        break;
    }
    N.pipesCount += 1, o("pipe count=%d opts=%j", N.pipesCount, t);
    var J = (!t || t.end !== !1) && r !== q.stdout && r !== q.stderr, k = J ? oe : de;
    N.endEmitted ? q.nextTick(k) : T.once("end", k), r.on("unpipe", V);
    function V(fe, ue) {
      o("onunpipe"), fe === T && ue && ue.hasUnpiped === !1 && (ue.hasUnpiped = !0, pr());
    }
    function oe() {
      o("onend"), r.end();
    }
    var A = te(T);
    r.on("drain", A);
    var pe = !1;
    function pr() {
      o("cleanup"), r.removeListener("close", we), r.removeListener("finish", Se), r.removeListener("drain", A), r.removeListener("error", me), r.removeListener("unpipe", V), T.removeListener("end", oe), T.removeListener("end", de), T.removeListener("data", yr), pe = !0, N.awaitDrain && (!r._writableState || r._writableState.needDrain) && A();
    }
    T.on("data", yr);
    function yr(fe) {
      o("ondata");
      var ue = r.write(fe);
      o("dest.write", ue), ue === !1 && ((N.pipesCount === 1 && N.pipes === r || N.pipesCount > 1 && B(N.pipes, r) !== -1) && !pe && (o("false write response, pause", N.awaitDrain), N.awaitDrain++), T.pause());
    }
    function me(fe) {
      o("onerror", fe), de(), r.removeListener("error", me), n(r, "error") === 0 && I(r, fe);
    }
    x(r, "error", me);
    function we() {
      r.removeListener("finish", Se), de();
    }
    r.once("close", we);
    function Se() {
      o("onfinish"), r.removeListener("close", we), de();
    }
    r.once("finish", Se);
    function de() {
      o("unpipe"), T.unpipe(r);
    }
    return r.emit("pipe", T), N.flowing || (o("pipe resume"), T.resume()), r;
  };
  function te(r) {
    return function() {
      var T = r._readableState;
      o("pipeOnDrain", T.awaitDrain), T.awaitDrain && T.awaitDrain--, T.awaitDrain === 0 && n(r, "data") && (T.flowing = !0, c(r));
    };
  }
  C.prototype.unpipe = function(r) {
    var t = this._readableState, T = {
      hasUnpiped: !1
    };
    if (t.pipesCount === 0)
      return this;
    if (t.pipesCount === 1)
      return r && r !== t.pipes ? this : (r || (r = t.pipes), t.pipes = null, t.pipesCount = 0, t.flowing = !1, r && r.emit("unpipe", this, T), this);
    if (!r) {
      var N = t.pipes, J = t.pipesCount;
      t.pipes = null, t.pipesCount = 0, t.flowing = !1;
      for (var k = 0; k < J; k++)
        N[k].emit("unpipe", this, {
          hasUnpiped: !1
        });
      return this;
    }
    var V = B(t.pipes, r);
    return V === -1 ? this : (t.pipes.splice(V, 1), t.pipesCount -= 1, t.pipesCount === 1 && (t.pipes = t.pipes[0]), r.emit("unpipe", this, T), this);
  }, C.prototype.on = function(r, t) {
    var T = i.prototype.on.call(this, r, t), N = this._readableState;
    return r === "data" ? (N.readableListening = this.listenerCount("readable") > 0, N.flowing !== !1 && this.resume()) : r === "readable" && !N.endEmitted && !N.readableListening && (N.readableListening = N.needReadable = !0, N.flowing = !1, N.emittedReadable = !1, o("on readable", N.length, N.reading), N.length ? Z(this) : N.reading || q.nextTick(m, this)), T;
  }, C.prototype.addListener = C.prototype.on, C.prototype.removeListener = function(r, t) {
    var T = i.prototype.removeListener.call(this, r, t);
    return r === "readable" && q.nextTick(p, this), T;
  }, C.prototype.removeAllListeners = function(r) {
    var t = i.prototype.removeAllListeners.apply(this, arguments);
    return (r === "readable" || r === void 0) && q.nextTick(p, this), t;
  };
  function p(r) {
    var t = r._readableState;
    t.readableListening = r.listenerCount("readable") > 0, t.resumeScheduled && !t.paused ? t.flowing = !0 : r.listenerCount("data") > 0 && r.resume();
  }
  function m(r) {
    o("readable nexttick read 0"), r.read(0);
  }
  C.prototype.resume = function() {
    var r = this._readableState;
    return r.flowing || (o("resume"), r.flowing = !r.readableListening, R(this, r)), r.paused = !1, this;
  };
  function R(r, t) {
    t.resumeScheduled || (t.resumeScheduled = !0, q.nextTick(F, r, t));
  }
  function F(r, t) {
    o("resume", t.reading), t.reading || r.read(0), t.resumeScheduled = !1, r.emit("resume"), c(r), t.flowing && !t.reading && r.read(0);
  }
  C.prototype.pause = function() {
    return o("call pause flowing=%j", this._readableState.flowing), this._readableState.flowing !== !1 && (o("pause"), this._readableState.flowing = !1, this.emit("pause")), this._readableState.paused = !0, this;
  };
  function c(r) {
    var t = r._readableState;
    for (o("flow", t.flowing); t.flowing && r.read() !== null; )
      ;
  }
  C.prototype.wrap = function(r) {
    var t = this, T = this._readableState, N = !1;
    r.on("end", function() {
      if (o("wrapped end"), T.decoder && !T.ended) {
        var V = T.decoder.end();
        V && V.length && t.push(V);
      }
      t.push(null);
    }), r.on("data", function(V) {
      if (o("wrapped data"), T.decoder && (V = T.decoder.write(V)), !(T.objectMode && V == null) && !(!T.objectMode && (!V || !V.length))) {
        var oe = t.push(V);
        oe || (N = !0, r.pause());
      }
    });
    for (var J in r)
      this[J] === void 0 && typeof r[J] == "function" && (this[J] = function(oe) {
        return function() {
          return r[oe].apply(r, arguments);
        };
      }(J));
    for (var k = 0; k < D.length; k++)
      r.on(D[k], this.emit.bind(this, D[k]));
    return this._read = function(V) {
      o("wrapped _read", V), N && (N = !1, r.resume());
    }, this;
  }, typeof Symbol == "function" && (C.prototype[Symbol.asyncIterator] = function() {
    return b === void 0 && (b = Gt()), b(this);
  }), Object.defineProperty(C.prototype, "readableHighWaterMark", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._readableState.highWaterMark;
    }
  }), Object.defineProperty(C.prototype, "readableBuffer", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._readableState && this._readableState.buffer;
    }
  }), Object.defineProperty(C.prototype, "readableFlowing", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._readableState.flowing;
    },
    set: function(t) {
      this._readableState && (this._readableState.flowing = t);
    }
  }), C._fromList = u, Object.defineProperty(C.prototype, "readableLength", {
    // making it explicit this property is not enumerable
    // because otherwise some prototype manipulation in
    // userland will fail
    enumerable: !1,
    get: function() {
      return this._readableState.length;
    }
  });
  function u(r, t) {
    if (t.length === 0)
      return null;
    var T;
    return t.objectMode ? T = t.buffer.shift() : !r || r >= t.length ? (t.decoder ? T = t.buffer.join("") : t.buffer.length === 1 ? T = t.buffer.first() : T = t.buffer.concat(t.length), t.buffer.clear()) : T = t.buffer.consume(r, t.decoder), T;
  }
  function w(r) {
    var t = r._readableState;
    o("endReadable", t.endEmitted), t.endEmitted || (t.ended = !0, q.nextTick(M, t, r));
  }
  function M(r, t) {
    if (o("endReadableNT", r.endEmitted, r.length), !r.endEmitted && r.length === 0 && (r.endEmitted = !0, t.readable = !1, t.emit("end"), r.autoDestroy)) {
      var T = t._writableState;
      (!T || T.autoDestroy && T.finished) && t.destroy();
    }
  }
  typeof Symbol == "function" && (C.from = function(r, t) {
    return P === void 0 && (P = Ht()), P(C, r, t);
  });
  function B(r, t) {
    for (var T = 0, N = r.length; T < N; T++)
      if (r[T] === t)
        return T;
    return -1;
  }
  return rr;
}
var tr, Xr;
function St() {
  if (Xr)
    return tr;
  Xr = 1, tr = l;
  var e = ce().codes, n = e.ERR_METHOD_NOT_IMPLEMENTED, i = e.ERR_MULTIPLE_CALLBACK, d = e.ERR_TRANSFORM_ALREADY_TRANSFORMING, y = e.ERR_TRANSFORM_WITH_LENGTH_0, s = le();
  ae(l, s);
  function f(v, S) {
    var O = this._transformState;
    O.transforming = !1;
    var L = O.writecb;
    if (L === null)
      return this.emit("error", new i());
    O.writechunk = null, O.writecb = null, S != null && this.push(S), L(v);
    var _ = this._readableState;
    _.reading = !1, (_.needReadable || _.length < _.highWaterMark) && this._read(_.highWaterMark);
  }
  function l(v) {
    if (!(this instanceof l))
      return new l(v);
    s.call(this, v), this._transformState = {
      afterTransform: f.bind(this),
      needTransform: !1,
      transforming: !1,
      writecb: null,
      writechunk: null,
      writeencoding: null
    }, this._readableState.needReadable = !0, this._readableState.sync = !1, v && (typeof v.transform == "function" && (this._transform = v.transform), typeof v.flush == "function" && (this._flush = v.flush)), this.on("prefinish", o);
  }
  function o() {
    var v = this;
    typeof this._flush == "function" && !this._readableState.destroyed ? this._flush(function(S, O) {
      h(v, S, O);
    }) : h(this, null, null);
  }
  l.prototype.push = function(v, S) {
    return this._transformState.needTransform = !1, s.prototype.push.call(this, v, S);
  }, l.prototype._transform = function(v, S, O) {
    O(new n("_transform()"));
  }, l.prototype._write = function(v, S, O) {
    var L = this._transformState;
    if (L.writecb = O, L.writechunk = v, L.writeencoding = S, !L.transforming) {
      var _ = this._readableState;
      (L.needTransform || _.needReadable || _.length < _.highWaterMark) && this._read(_.highWaterMark);
    }
  }, l.prototype._read = function(v) {
    var S = this._transformState;
    S.writechunk !== null && !S.transforming ? (S.transforming = !0, this._transform(S.writechunk, S.writeencoding, S.afterTransform)) : S.needTransform = !0;
  }, l.prototype._destroy = function(v, S) {
    s.prototype._destroy.call(this, v, function(O) {
      S(O);
    });
  };
  function h(v, S, O) {
    if (S)
      return v.emit("error", S);
    if (O != null && v.push(O), v._writableState.length)
      throw new y();
    if (v._transformState.transforming)
      throw new d();
    return v.push(null);
  }
  return tr;
}
var nr, et;
function $t() {
  if (et)
    return nr;
  et = 1, nr = n;
  var e = St();
  ae(n, e);
  function n(i) {
    if (!(this instanceof n))
      return new n(i);
    e.call(this, i);
  }
  return n.prototype._transform = function(i, d, y) {
    y(null, i);
  }, nr;
}
var ir, rt;
function Vt() {
  if (rt)
    return ir;
  rt = 1;
  var e;
  function n(O) {
    var L = !1;
    return function() {
      L || (L = !0, O.apply(void 0, arguments));
    };
  }
  var i = ce().codes, d = i.ERR_MISSING_ARGS, y = i.ERR_STREAM_DESTROYED;
  function s(O) {
    if (O)
      throw O;
  }
  function f(O) {
    return O.setHeader && typeof O.abort == "function";
  }
  function l(O, L, _, j) {
    j = n(j);
    var E = !1;
    O.on("close", function() {
      E = !0;
    }), e === void 0 && (e = cr()), e(O, {
      readable: L,
      writable: _
    }, function(g) {
      if (g)
        return j(g);
      E = !0, j();
    });
    var a = !1;
    return function(g) {
      if (!E && !a) {
        if (a = !0, f(O))
          return O.abort();
        if (typeof O.destroy == "function")
          return O.destroy();
        j(g || new y("pipe"));
      }
    };
  }
  function o(O) {
    O();
  }
  function h(O, L) {
    return O.pipe(L);
  }
  function v(O) {
    return !O.length || typeof O[O.length - 1] != "function" ? s : O.pop();
  }
  function S() {
    for (var O = arguments.length, L = new Array(O), _ = 0; _ < O; _++)
      L[_] = arguments[_];
    var j = v(L);
    if (Array.isArray(L[0]) && (L = L[0]), L.length < 2)
      throw new d("streams");
    var E, a = L.map(function(g, b) {
      var P = b < L.length - 1, I = b > 0;
      return l(g, P, I, function(D) {
        E || (E = D), D && a.forEach(o), !P && (a.forEach(o), j(E));
      });
    });
    return L.reduce(h);
  }
  return ir = S, ir;
}
var zt = ne, dr = fr.EventEmitter, Kt = ae;
Kt(ne, dr);
ne.Readable = wt();
ne.Writable = mt();
ne.Duplex = le();
ne.Transform = St();
ne.PassThrough = $t();
ne.finished = cr();
ne.pipeline = Vt();
ne.Stream = ne;
function ne() {
  dr.call(this);
}
ne.prototype.pipe = function(e, n) {
  var i = this;
  function d(v) {
    e.writable && e.write(v) === !1 && i.pause && i.pause();
  }
  i.on("data", d);
  function y() {
    i.readable && i.resume && i.resume();
  }
  e.on("drain", y), !e._isStdio && (!n || n.end !== !1) && (i.on("end", f), i.on("close", l));
  var s = !1;
  function f() {
    s || (s = !0, e.end());
  }
  function l() {
    s || (s = !0, typeof e.destroy == "function" && e.destroy());
  }
  function o(v) {
    if (h(), dr.listenerCount(this, "error") === 0)
      throw v;
  }
  i.on("error", o), e.on("error", o);
  function h() {
    i.removeListener("data", d), e.removeListener("drain", y), i.removeListener("end", f), i.removeListener("close", l), i.removeListener("error", o), e.removeListener("error", o), i.removeListener("end", h), i.removeListener("close", h), e.removeListener("close", h);
  }
  return i.on("end", h), i.on("close", h), e.on("close", h), e.emit("pipe", i), e;
};
const Yt = /* @__PURE__ */ _t(zt);
export {
  Yt as S,
  wt as a,
  mt as b,
  le as c,
  St as d,
  $t as e,
  cr as f,
  Vt as g,
  gt as h,
  ae as i,
  sr as j,
  Nt as k,
  Kr as r,
  zt as s
};
