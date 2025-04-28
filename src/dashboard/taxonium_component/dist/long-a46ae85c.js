import { Q as D } from "./index-8ee3baa4.js";
var H = r, v = null;
try {
  v = new WebAssembly.Instance(new WebAssembly.Module(new Uint8Array([
    0,
    97,
    115,
    109,
    1,
    0,
    0,
    0,
    1,
    13,
    2,
    96,
    0,
    1,
    127,
    96,
    4,
    127,
    127,
    127,
    127,
    1,
    127,
    3,
    7,
    6,
    0,
    1,
    1,
    1,
    1,
    1,
    6,
    6,
    1,
    127,
    1,
    65,
    0,
    11,
    7,
    50,
    6,
    3,
    109,
    117,
    108,
    0,
    1,
    5,
    100,
    105,
    118,
    95,
    115,
    0,
    2,
    5,
    100,
    105,
    118,
    95,
    117,
    0,
    3,
    5,
    114,
    101,
    109,
    95,
    115,
    0,
    4,
    5,
    114,
    101,
    109,
    95,
    117,
    0,
    5,
    8,
    103,
    101,
    116,
    95,
    104,
    105,
    103,
    104,
    0,
    0,
    10,
    191,
    1,
    6,
    4,
    0,
    35,
    0,
    11,
    36,
    1,
    1,
    126,
    32,
    0,
    173,
    32,
    1,
    173,
    66,
    32,
    134,
    132,
    32,
    2,
    173,
    32,
    3,
    173,
    66,
    32,
    134,
    132,
    126,
    34,
    4,
    66,
    32,
    135,
    167,
    36,
    0,
    32,
    4,
    167,
    11,
    36,
    1,
    1,
    126,
    32,
    0,
    173,
    32,
    1,
    173,
    66,
    32,
    134,
    132,
    32,
    2,
    173,
    32,
    3,
    173,
    66,
    32,
    134,
    132,
    127,
    34,
    4,
    66,
    32,
    135,
    167,
    36,
    0,
    32,
    4,
    167,
    11,
    36,
    1,
    1,
    126,
    32,
    0,
    173,
    32,
    1,
    173,
    66,
    32,
    134,
    132,
    32,
    2,
    173,
    32,
    3,
    173,
    66,
    32,
    134,
    132,
    128,
    34,
    4,
    66,
    32,
    135,
    167,
    36,
    0,
    32,
    4,
    167,
    11,
    36,
    1,
    1,
    126,
    32,
    0,
    173,
    32,
    1,
    173,
    66,
    32,
    134,
    132,
    32,
    2,
    173,
    32,
    3,
    173,
    66,
    32,
    134,
    132,
    129,
    34,
    4,
    66,
    32,
    135,
    167,
    36,
    0,
    32,
    4,
    167,
    11,
    36,
    1,
    1,
    126,
    32,
    0,
    173,
    32,
    1,
    173,
    66,
    32,
    134,
    132,
    32,
    2,
    173,
    32,
    3,
    173,
    66,
    32,
    134,
    132,
    130,
    34,
    4,
    66,
    32,
    135,
    167,
    36,
    0,
    32,
    4,
    167,
    11
  ])), {}).exports;
} catch {
}
function r(i, t, e) {
  this.low = i | 0, this.high = t | 0, this.unsigned = !!e;
}
r.prototype.__isLong__;
Object.defineProperty(r.prototype, "__isLong__", { value: !0 });
function w(i) {
  return (i && i.__isLong__) === !0;
}
r.isLong = w;
var A = {}, B = {};
function U(i, t) {
  var e, h, g;
  return t ? (i >>>= 0, (g = 0 <= i && i < 256) && (h = B[i], h) ? h : (e = s(i, (i | 0) < 0 ? -1 : 0, !0), g && (B[i] = e), e)) : (i |= 0, (g = -128 <= i && i < 128) && (h = A[i], h) ? h : (e = s(i, i < 0 ? -1 : 0, !1), g && (A[i] = e), e));
}
r.fromInt = U;
function L(i, t) {
  if (isNaN(i))
    return t ? m : q;
  if (t) {
    if (i < 0)
      return m;
    if (i >= P)
      return V;
  } else {
    if (i <= -W)
      return c;
    if (i + 1 >= W)
      return R;
  }
  return i < 0 ? L(-i, t).neg() : s(i % x | 0, i / x | 0, t);
}
r.fromNumber = L;
function s(i, t, e) {
  return new r(i, t, e);
}
r.fromBits = s;
var b = Math.pow;
function d(i, t, e) {
  if (i.length === 0)
    throw Error("empty string");
  if (i === "NaN" || i === "Infinity" || i === "+Infinity" || i === "-Infinity")
    return q;
  if (typeof t == "number" ? (e = t, t = !1) : t = !!t, e = e || 10, e < 2 || 36 < e)
    throw RangeError("radix");
  var h;
  if ((h = i.indexOf("-")) > 0)
    throw Error("interior hyphen");
  if (h === 0)
    return d(i.substring(1), t, e).neg();
  for (var g = L(b(e, 8)), f = q, o = 0; o < i.length; o += 8) {
    var l = Math.min(8, i.length - o), _ = parseInt(i.substring(o, o + l), e);
    if (l < 8) {
      var F = L(b(e, l));
      f = f.mul(F).add(L(_));
    } else
      f = f.mul(g), f = f.add(L(_));
  }
  return f.unsigned = t, f;
}
r.fromString = d;
function O(i, t) {
  return typeof i == "number" ? L(i, t) : typeof i == "string" ? d(i, t) : s(i.low, i.high, typeof t == "boolean" ? t : i.unsigned);
}
r.fromValue = O;
var M = 65536, C = 1 << 24, x = M * M, P = x * x, W = P / 2, Z = U(C), q = U(0);
r.ZERO = q;
var m = U(0, !0);
r.UZERO = m;
var T = U(1);
r.ONE = T;
var S = U(1, !0);
r.UONE = S;
var I = U(-1);
r.NEG_ONE = I;
var R = s(-1, 2147483647, !1);
r.MAX_VALUE = R;
var V = s(-1, -1, !0);
r.MAX_UNSIGNED_VALUE = V;
var c = s(0, -2147483648, !1);
r.MIN_VALUE = c;
var n = r.prototype;
n.toInt = function() {
  return this.unsigned ? this.low >>> 0 : this.low;
};
n.toNumber = function() {
  return this.unsigned ? (this.high >>> 0) * x + (this.low >>> 0) : this.high * x + (this.low >>> 0);
};
n.toString = function(t) {
  if (t = t || 10, t < 2 || 36 < t)
    throw RangeError("radix");
  if (this.isZero())
    return "0";
  if (this.isNegative())
    if (this.eq(c)) {
      var e = L(t), h = this.div(e), g = h.mul(e).sub(this);
      return h.toString(t) + g.toInt().toString(t);
    } else
      return "-" + this.neg().toString(t);
  for (var f = L(b(t, 6), this.unsigned), o = this, l = ""; ; ) {
    var _ = o.div(f), F = o.sub(_.mul(f)).toInt() >>> 0, u = F.toString(t);
    if (o = _, o.isZero())
      return u + l;
    for (; u.length < 6; )
      u = "0" + u;
    l = "" + u + l;
  }
};
n.getHighBits = function() {
  return this.high;
};
n.getHighBitsUnsigned = function() {
  return this.high >>> 0;
};
n.getLowBits = function() {
  return this.low;
};
n.getLowBitsUnsigned = function() {
  return this.low >>> 0;
};
n.getNumBitsAbs = function() {
  if (this.isNegative())
    return this.eq(c) ? 64 : this.neg().getNumBitsAbs();
  for (var t = this.high != 0 ? this.high : this.low, e = 31; e > 0 && !(t & 1 << e); e--)
    ;
  return this.high != 0 ? e + 33 : e + 1;
};
n.isZero = function() {
  return this.high === 0 && this.low === 0;
};
n.eqz = n.isZero;
n.isNegative = function() {
  return !this.unsigned && this.high < 0;
};
n.isPositive = function() {
  return this.unsigned || this.high >= 0;
};
n.isOdd = function() {
  return (this.low & 1) === 1;
};
n.isEven = function() {
  return (this.low & 1) === 0;
};
n.equals = function(t) {
  return w(t) || (t = O(t)), this.unsigned !== t.unsigned && this.high >>> 31 === 1 && t.high >>> 31 === 1 ? !1 : this.high === t.high && this.low === t.low;
};
n.eq = n.equals;
n.notEquals = function(t) {
  return !this.eq(
    /* validates */
    t
  );
};
n.neq = n.notEquals;
n.ne = n.notEquals;
n.lessThan = function(t) {
  return this.comp(
    /* validates */
    t
  ) < 0;
};
n.lt = n.lessThan;
n.lessThanOrEqual = function(t) {
  return this.comp(
    /* validates */
    t
  ) <= 0;
};
n.lte = n.lessThanOrEqual;
n.le = n.lessThanOrEqual;
n.greaterThan = function(t) {
  return this.comp(
    /* validates */
    t
  ) > 0;
};
n.gt = n.greaterThan;
n.greaterThanOrEqual = function(t) {
  return this.comp(
    /* validates */
    t
  ) >= 0;
};
n.gte = n.greaterThanOrEqual;
n.ge = n.greaterThanOrEqual;
n.compare = function(t) {
  if (w(t) || (t = O(t)), this.eq(t))
    return 0;
  var e = this.isNegative(), h = t.isNegative();
  return e && !h ? -1 : !e && h ? 1 : this.unsigned ? t.high >>> 0 > this.high >>> 0 || t.high === this.high && t.low >>> 0 > this.low >>> 0 ? -1 : 1 : this.sub(t).isNegative() ? -1 : 1;
};
n.comp = n.compare;
n.negate = function() {
  return !this.unsigned && this.eq(c) ? c : this.not().add(T);
};
n.neg = n.negate;
n.add = function(t) {
  w(t) || (t = O(t));
  var e = this.high >>> 16, h = this.high & 65535, g = this.low >>> 16, f = this.low & 65535, o = t.high >>> 16, l = t.high & 65535, _ = t.low >>> 16, F = t.low & 65535, u = 0, N = 0, a = 0, E = 0;
  return E += f + F, a += E >>> 16, E &= 65535, a += g + _, N += a >>> 16, a &= 65535, N += h + l, u += N >>> 16, N &= 65535, u += e + o, u &= 65535, s(a << 16 | E, u << 16 | N, this.unsigned);
};
n.subtract = function(t) {
  return w(t) || (t = O(t)), this.add(t.neg());
};
n.sub = n.subtract;
n.multiply = function(t) {
  if (this.isZero())
    return q;
  if (w(t) || (t = O(t)), v) {
    var e = v.mul(
      this.low,
      this.high,
      t.low,
      t.high
    );
    return s(e, v.get_high(), this.unsigned);
  }
  if (t.isZero())
    return q;
  if (this.eq(c))
    return t.isOdd() ? c : q;
  if (t.eq(c))
    return this.isOdd() ? c : q;
  if (this.isNegative())
    return t.isNegative() ? this.neg().mul(t.neg()) : this.neg().mul(t).neg();
  if (t.isNegative())
    return this.mul(t.neg()).neg();
  if (this.lt(Z) && t.lt(Z))
    return L(this.toNumber() * t.toNumber(), this.unsigned);
  var h = this.high >>> 16, g = this.high & 65535, f = this.low >>> 16, o = this.low & 65535, l = t.high >>> 16, _ = t.high & 65535, F = t.low >>> 16, u = t.low & 65535, N = 0, a = 0, E = 0, y = 0;
  return y += o * u, E += y >>> 16, y &= 65535, E += f * u, a += E >>> 16, E &= 65535, E += o * F, a += E >>> 16, E &= 65535, a += g * u, N += a >>> 16, a &= 65535, a += f * F, N += a >>> 16, a &= 65535, a += o * _, N += a >>> 16, a &= 65535, N += h * u + g * F + f * _ + o * l, N &= 65535, s(E << 16 | y, N << 16 | a, this.unsigned);
};
n.mul = n.multiply;
n.divide = function(t) {
  if (w(t) || (t = O(t)), t.isZero())
    throw Error("division by zero");
  if (v) {
    if (!this.unsigned && this.high === -2147483648 && t.low === -1 && t.high === -1)
      return this;
    var e = (this.unsigned ? v.div_u : v.div_s)(
      this.low,
      this.high,
      t.low,
      t.high
    );
    return s(e, v.get_high(), this.unsigned);
  }
  if (this.isZero())
    return this.unsigned ? m : q;
  var h, g, f;
  if (this.unsigned) {
    if (t.unsigned || (t = t.toUnsigned()), t.gt(this))
      return m;
    if (t.gt(this.shru(1)))
      return S;
    f = m;
  } else {
    if (this.eq(c)) {
      if (t.eq(T) || t.eq(I))
        return c;
      if (t.eq(c))
        return T;
      var o = this.shr(1);
      return h = o.div(t).shl(1), h.eq(q) ? t.isNegative() ? T : I : (g = this.sub(t.mul(h)), f = h.add(g.div(t)), f);
    } else if (t.eq(c))
      return this.unsigned ? m : q;
    if (this.isNegative())
      return t.isNegative() ? this.neg().div(t.neg()) : this.neg().div(t).neg();
    if (t.isNegative())
      return this.div(t.neg()).neg();
    f = q;
  }
  for (g = this; g.gte(t); ) {
    h = Math.max(1, Math.floor(g.toNumber() / t.toNumber()));
    for (var l = Math.ceil(Math.log(h) / Math.LN2), _ = l <= 48 ? 1 : b(2, l - 48), F = L(h), u = F.mul(t); u.isNegative() || u.gt(g); )
      h -= _, F = L(h, this.unsigned), u = F.mul(t);
    F.isZero() && (F = T), f = f.add(F), g = g.sub(u);
  }
  return f;
};
n.div = n.divide;
n.modulo = function(t) {
  if (w(t) || (t = O(t)), v) {
    var e = (this.unsigned ? v.rem_u : v.rem_s)(
      this.low,
      this.high,
      t.low,
      t.high
    );
    return s(e, v.get_high(), this.unsigned);
  }
  return this.sub(this.div(t).mul(t));
};
n.mod = n.modulo;
n.rem = n.modulo;
n.not = function() {
  return s(~this.low, ~this.high, this.unsigned);
};
n.and = function(t) {
  return w(t) || (t = O(t)), s(this.low & t.low, this.high & t.high, this.unsigned);
};
n.or = function(t) {
  return w(t) || (t = O(t)), s(this.low | t.low, this.high | t.high, this.unsigned);
};
n.xor = function(t) {
  return w(t) || (t = O(t)), s(this.low ^ t.low, this.high ^ t.high, this.unsigned);
};
n.shiftLeft = function(t) {
  return w(t) && (t = t.toInt()), (t &= 63) === 0 ? this : t < 32 ? s(this.low << t, this.high << t | this.low >>> 32 - t, this.unsigned) : s(0, this.low << t - 32, this.unsigned);
};
n.shl = n.shiftLeft;
n.shiftRight = function(t) {
  return w(t) && (t = t.toInt()), (t &= 63) === 0 ? this : t < 32 ? s(this.low >>> t | this.high << 32 - t, this.high >> t, this.unsigned) : s(this.high >> t - 32, this.high >= 0 ? 0 : -1, this.unsigned);
};
n.shr = n.shiftRight;
n.shiftRightUnsigned = function(t) {
  if (w(t) && (t = t.toInt()), t &= 63, t === 0)
    return this;
  var e = this.high;
  if (t < 32) {
    var h = this.low;
    return s(h >>> t | e << 32 - t, e >>> t, this.unsigned);
  } else
    return t === 32 ? s(e, 0, this.unsigned) : s(e >>> t - 32, 0, this.unsigned);
};
n.shru = n.shiftRightUnsigned;
n.shr_u = n.shiftRightUnsigned;
n.toSigned = function() {
  return this.unsigned ? s(this.low, this.high, !1) : this;
};
n.toUnsigned = function() {
  return this.unsigned ? this : s(this.low, this.high, !0);
};
n.toBytes = function(t) {
  return t ? this.toBytesLE() : this.toBytesBE();
};
n.toBytesLE = function() {
  var t = this.high, e = this.low;
  return [
    e & 255,
    e >>> 8 & 255,
    e >>> 16 & 255,
    e >>> 24,
    t & 255,
    t >>> 8 & 255,
    t >>> 16 & 255,
    t >>> 24
  ];
};
n.toBytesBE = function() {
  var t = this.high, e = this.low;
  return [
    t >>> 24,
    t >>> 16 & 255,
    t >>> 8 & 255,
    t & 255,
    e >>> 24,
    e >>> 16 & 255,
    e >>> 8 & 255,
    e & 255
  ];
};
r.fromBytes = function(t, e, h) {
  return h ? r.fromBytesLE(t, e) : r.fromBytesBE(t, e);
};
r.fromBytesLE = function(t, e) {
  return new r(
    t[0] | t[1] << 8 | t[2] << 16 | t[3] << 24,
    t[4] | t[5] << 8 | t[6] << 16 | t[7] << 24,
    e
  );
};
r.fromBytesBE = function(t, e) {
  return new r(
    t[4] << 24 | t[5] << 16 | t[6] << 8 | t[7],
    t[0] << 24 | t[1] << 16 | t[2] << 8 | t[3],
    e
  );
};
const X = /* @__PURE__ */ D(H);
export {
  X as L,
  H as l
};
