import p, { useState as ee } from "react";
import { K as yt, N as bt, O as Fe, Q as St, o as Dt, q as He, D as _t, S as wt, a as Ot, I as Pt, U as Et, J as xt, P as Ct, r as De, i as _e, j as we, V as Mt, t as q, b as Tt, B as Y, c as Nt, W as Rt, X as jt, Y as At } from "./index-8ee3baa4.js";
import qe from "react-dom";
import { a as kt, b as $t } from "./ColorPicker-affa19c8.js";
import "./index-090e001f.js";
var ce = { exports: {} }, Ve = {};
const It = /* @__PURE__ */ yt(bt);
var _ = {}, W = {};
Object.defineProperty(W, "__esModule", {
  value: !0
});
W.dontSetMe = Bt;
W.findInArray = Wt;
W.int = Yt;
W.isFunction = Ut;
W.isNum = Lt;
function Wt(e, t) {
  for (var n = 0, r = e.length; n < r; n++)
    if (t.apply(t, [e[n], n, e]))
      return e[n];
}
function Ut(e) {
  return typeof e == "function" || Object.prototype.toString.call(e) === "[object Function]";
}
function Lt(e) {
  return typeof e == "number" && !isNaN(e);
}
function Yt(e) {
  return parseInt(e, 10);
}
function Bt(e, t, n) {
  if (e[t])
    return new Error("Invalid prop ".concat(t, " passed to ").concat(n, " - do not set this, set it on the child."));
}
var V = {};
Object.defineProperty(V, "__esModule", {
  value: !0
});
V.browserPrefixToKey = Ke;
V.browserPrefixToStyle = Xt;
V.default = void 0;
V.getPrefix = Ge;
var he = ["Moz", "Webkit", "O", "ms"];
function Ge() {
  var e, t, n = arguments.length > 0 && arguments[0] !== void 0 ? arguments[0] : "transform";
  if (typeof window > "u")
    return "";
  var r = (e = window.document) === null || e === void 0 || (t = e.documentElement) === null || t === void 0 ? void 0 : t.style;
  if (!r || n in r)
    return "";
  for (var o = 0; o < he.length; o++)
    if (Ke(n, he[o]) in r)
      return he[o];
  return "";
}
function Ke(e, t) {
  return t ? "".concat(t).concat(Ft(e)) : e;
}
function Xt(e, t) {
  return t ? "-".concat(t.toLowerCase(), "-").concat(e) : e;
}
function Ft(e) {
  for (var t = "", n = !0, r = 0; r < e.length; r++)
    n ? (t += e[r].toUpperCase(), n = !1) : e[r] === "-" ? n = !0 : t += e[r];
  return t;
}
var Ht = Ge();
V.default = Ht;
function be(e) {
  "@babel/helpers - typeof";
  return be = typeof Symbol == "function" && typeof Symbol.iterator == "symbol" ? function(t) {
    return typeof t;
  } : function(t) {
    return t && typeof Symbol == "function" && t.constructor === Symbol && t !== Symbol.prototype ? "symbol" : typeof t;
  }, be(e);
}
Object.defineProperty(_, "__esModule", {
  value: !0
});
_.addClassName = et;
_.addEvent = Gt;
_.addUserSelectStyles = or;
_.createCSSTransform = tr;
_.createSVGTransform = rr;
_.getTouch = nr;
_.getTouchIdentifier = ar;
_.getTranslation = Oe;
_.innerHeight = Qt;
_.innerWidth = Zt;
_.matchesSelector = Ze;
_.matchesSelectorAndParentsTo = Vt;
_.offsetXYFromParent = er;
_.outerHeight = zt;
_.outerWidth = Jt;
_.removeClassName = tt;
_.removeEvent = Kt;
_.removeUserSelectStyles = ir;
var R = W, ke = qt(V);
function ze(e) {
  if (typeof WeakMap != "function")
    return null;
  var t = /* @__PURE__ */ new WeakMap(), n = /* @__PURE__ */ new WeakMap();
  return (ze = function(o) {
    return o ? n : t;
  })(e);
}
function qt(e, t) {
  if (!t && e && e.__esModule)
    return e;
  if (e === null || be(e) !== "object" && typeof e != "function")
    return { default: e };
  var n = ze(t);
  if (n && n.has(e))
    return n.get(e);
  var r = {}, o = Object.defineProperty && Object.getOwnPropertyDescriptor;
  for (var i in e)
    if (i !== "default" && Object.prototype.hasOwnProperty.call(e, i)) {
      var f = o ? Object.getOwnPropertyDescriptor(e, i) : null;
      f && (f.get || f.set) ? Object.defineProperty(r, i, f) : r[i] = e[i];
    }
  return r.default = e, n && n.set(e, r), r;
}
function $e(e, t) {
  var n = Object.keys(e);
  if (Object.getOwnPropertySymbols) {
    var r = Object.getOwnPropertySymbols(e);
    t && (r = r.filter(function(o) {
      return Object.getOwnPropertyDescriptor(e, o).enumerable;
    })), n.push.apply(n, r);
  }
  return n;
}
function Je(e) {
  for (var t = 1; t < arguments.length; t++) {
    var n = arguments[t] != null ? arguments[t] : {};
    t % 2 ? $e(Object(n), !0).forEach(function(r) {
      Qe(e, r, n[r]);
    }) : Object.getOwnPropertyDescriptors ? Object.defineProperties(e, Object.getOwnPropertyDescriptors(n)) : $e(Object(n)).forEach(function(r) {
      Object.defineProperty(e, r, Object.getOwnPropertyDescriptor(n, r));
    });
  }
  return e;
}
function Qe(e, t, n) {
  return t in e ? Object.defineProperty(e, t, { value: n, enumerable: !0, configurable: !0, writable: !0 }) : e[t] = n, e;
}
var le = "";
function Ze(e, t) {
  return le || (le = (0, R.findInArray)(["matches", "webkitMatchesSelector", "mozMatchesSelector", "msMatchesSelector", "oMatchesSelector"], function(n) {
    return (0, R.isFunction)(e[n]);
  })), (0, R.isFunction)(e[le]) ? e[le](t) : !1;
}
function Vt(e, t, n) {
  var r = e;
  do {
    if (Ze(r, t))
      return !0;
    if (r === n)
      return !1;
    r = r.parentNode;
  } while (r);
  return !1;
}
function Gt(e, t, n, r) {
  if (e) {
    var o = Je({
      capture: !0
    }, r);
    e.addEventListener ? e.addEventListener(t, n, o) : e.attachEvent ? e.attachEvent("on" + t, n) : e["on" + t] = n;
  }
}
function Kt(e, t, n, r) {
  if (e) {
    var o = Je({
      capture: !0
    }, r);
    e.removeEventListener ? e.removeEventListener(t, n, o) : e.detachEvent ? e.detachEvent("on" + t, n) : e["on" + t] = null;
  }
}
function zt(e) {
  var t = e.clientHeight, n = e.ownerDocument.defaultView.getComputedStyle(e);
  return t += (0, R.int)(n.borderTopWidth), t += (0, R.int)(n.borderBottomWidth), t;
}
function Jt(e) {
  var t = e.clientWidth, n = e.ownerDocument.defaultView.getComputedStyle(e);
  return t += (0, R.int)(n.borderLeftWidth), t += (0, R.int)(n.borderRightWidth), t;
}
function Qt(e) {
  var t = e.clientHeight, n = e.ownerDocument.defaultView.getComputedStyle(e);
  return t -= (0, R.int)(n.paddingTop), t -= (0, R.int)(n.paddingBottom), t;
}
function Zt(e) {
  var t = e.clientWidth, n = e.ownerDocument.defaultView.getComputedStyle(e);
  return t -= (0, R.int)(n.paddingLeft), t -= (0, R.int)(n.paddingRight), t;
}
function er(e, t, n) {
  var r = t === t.ownerDocument.body, o = r ? {
    left: 0,
    top: 0
  } : t.getBoundingClientRect(), i = (e.clientX + t.scrollLeft - o.left) / n, f = (e.clientY + t.scrollTop - o.top) / n;
  return {
    x: i,
    y: f
  };
}
function tr(e, t) {
  var n = Oe(e, t, "px");
  return Qe({}, (0, ke.browserPrefixToKey)("transform", ke.default), n);
}
function rr(e, t) {
  var n = Oe(e, t, "");
  return n;
}
function Oe(e, t, n) {
  var r = e.x, o = e.y, i = "translate(".concat(r).concat(n, ",").concat(o).concat(n, ")");
  if (t) {
    var f = "".concat(typeof t.x == "string" ? t.x : t.x + n), d = "".concat(typeof t.y == "string" ? t.y : t.y + n);
    i = "translate(".concat(f, ", ").concat(d, ")") + i;
  }
  return i;
}
function nr(e, t) {
  return e.targetTouches && (0, R.findInArray)(e.targetTouches, function(n) {
    return t === n.identifier;
  }) || e.changedTouches && (0, R.findInArray)(e.changedTouches, function(n) {
    return t === n.identifier;
  });
}
function ar(e) {
  if (e.targetTouches && e.targetTouches[0])
    return e.targetTouches[0].identifier;
  if (e.changedTouches && e.changedTouches[0])
    return e.changedTouches[0].identifier;
}
function or(e) {
  if (e) {
    var t = e.getElementById("react-draggable-style-el");
    t || (t = e.createElement("style"), t.type = "text/css", t.id = "react-draggable-style-el", t.innerHTML = `.react-draggable-transparent-selection *::-moz-selection {all: inherit;}
`, t.innerHTML += `.react-draggable-transparent-selection *::selection {all: inherit;}
`, e.getElementsByTagName("head")[0].appendChild(t)), e.body && et(e.body, "react-draggable-transparent-selection");
  }
}
function ir(e) {
  if (e)
    try {
      if (e.body && tt(e.body, "react-draggable-transparent-selection"), e.selection)
        e.selection.empty();
      else {
        var t = (e.defaultView || window).getSelection();
        t && t.type !== "Caret" && t.removeAllRanges();
      }
    } catch {
    }
}
function et(e, t) {
  e.classList ? e.classList.add(t) : e.className.match(new RegExp("(?:^|\\s)".concat(t, "(?!\\S)"))) || (e.className += " ".concat(t));
}
function tt(e, t) {
  e.classList ? e.classList.remove(t) : e.className = e.className.replace(new RegExp("(?:^|\\s)".concat(t, "(?!\\S)"), "g"), "");
}
var U = {};
Object.defineProperty(U, "__esModule", {
  value: !0
});
U.canDragX = sr;
U.canDragY = cr;
U.createCoreData = dr;
U.createDraggableData = pr;
U.getBoundPosition = lr;
U.getControlPosition = fr;
U.snapToGrid = ur;
var N = W, z = _;
function lr(e, t, n) {
  if (!e.props.bounds)
    return [t, n];
  var r = e.props.bounds;
  r = typeof r == "string" ? r : gr(r);
  var o = Pe(e);
  if (typeof r == "string") {
    var i = o.ownerDocument, f = i.defaultView, d;
    if (r === "parent" ? d = o.parentNode : d = i.querySelector(r), !(d instanceof f.HTMLElement))
      throw new Error('Bounds selector "' + r + '" could not find an element.');
    var v = d, b = f.getComputedStyle(o), S = f.getComputedStyle(v);
    r = {
      left: -o.offsetLeft + (0, N.int)(S.paddingLeft) + (0, N.int)(b.marginLeft),
      top: -o.offsetTop + (0, N.int)(S.paddingTop) + (0, N.int)(b.marginTop),
      right: (0, z.innerWidth)(v) - (0, z.outerWidth)(o) - o.offsetLeft + (0, N.int)(S.paddingRight) - (0, N.int)(b.marginRight),
      bottom: (0, z.innerHeight)(v) - (0, z.outerHeight)(o) - o.offsetTop + (0, N.int)(S.paddingBottom) - (0, N.int)(b.marginBottom)
    };
  }
  return (0, N.isNum)(r.right) && (t = Math.min(t, r.right)), (0, N.isNum)(r.bottom) && (n = Math.min(n, r.bottom)), (0, N.isNum)(r.left) && (t = Math.max(t, r.left)), (0, N.isNum)(r.top) && (n = Math.max(n, r.top)), [t, n];
}
function ur(e, t, n) {
  var r = Math.round(t / e[0]) * e[0], o = Math.round(n / e[1]) * e[1];
  return [r, o];
}
function sr(e) {
  return e.props.axis === "both" || e.props.axis === "x";
}
function cr(e) {
  return e.props.axis === "both" || e.props.axis === "y";
}
function fr(e, t, n) {
  var r = typeof t == "number" ? (0, z.getTouch)(e, t) : null;
  if (typeof t == "number" && !r)
    return null;
  var o = Pe(n), i = n.props.offsetParent || o.offsetParent || o.ownerDocument.body;
  return (0, z.offsetXYFromParent)(r || e, i, n.props.scale);
}
function dr(e, t, n) {
  var r = e.state, o = !(0, N.isNum)(r.lastX), i = Pe(e);
  return o ? {
    node: i,
    deltaX: 0,
    deltaY: 0,
    lastX: t,
    lastY: n,
    x: t,
    y: n
  } : {
    node: i,
    deltaX: t - r.lastX,
    deltaY: n - r.lastY,
    lastX: r.lastX,
    lastY: r.lastY,
    x: t,
    y: n
  };
}
function pr(e, t) {
  var n = e.props.scale;
  return {
    node: t.node,
    x: e.state.x + t.deltaX / n,
    y: e.state.y + t.deltaY / n,
    deltaX: t.deltaX / n,
    deltaY: t.deltaY / n,
    lastX: e.state.x,
    lastY: e.state.y
  };
}
function gr(e) {
  return {
    left: e.left,
    top: e.top,
    right: e.right,
    bottom: e.bottom
  };
}
function Pe(e) {
  var t = e.findDOMNode();
  if (!t)
    throw new Error("<DraggableCore>: Unmounted during event!");
  return t;
}
var fe = {}, de = {};
Object.defineProperty(de, "__esModule", {
  value: !0
});
de.default = mr;
function mr() {
}
function ue(e) {
  "@babel/helpers - typeof";
  return ue = typeof Symbol == "function" && typeof Symbol.iterator == "symbol" ? function(t) {
    return typeof t;
  } : function(t) {
    return t && typeof Symbol == "function" && t.constructor === Symbol && t !== Symbol.prototype ? "symbol" : typeof t;
  }, ue(e);
}
Object.defineProperty(fe, "__esModule", {
  value: !0
});
fe.default = void 0;
var ve = vr(p), T = Ee(Fe), hr = Ee(qe), x = _, X = U, ye = W, Z = Ee(de);
function Ee(e) {
  return e && e.__esModule ? e : { default: e };
}
function rt(e) {
  if (typeof WeakMap != "function")
    return null;
  var t = /* @__PURE__ */ new WeakMap(), n = /* @__PURE__ */ new WeakMap();
  return (rt = function(o) {
    return o ? n : t;
  })(e);
}
function vr(e, t) {
  if (!t && e && e.__esModule)
    return e;
  if (e === null || ue(e) !== "object" && typeof e != "function")
    return { default: e };
  var n = rt(t);
  if (n && n.has(e))
    return n.get(e);
  var r = {}, o = Object.defineProperty && Object.getOwnPropertyDescriptor;
  for (var i in e)
    if (i !== "default" && Object.prototype.hasOwnProperty.call(e, i)) {
      var f = o ? Object.getOwnPropertyDescriptor(e, i) : null;
      f && (f.get || f.set) ? Object.defineProperty(r, i, f) : r[i] = e[i];
    }
  return r.default = e, n && n.set(e, r), r;
}
function Ie(e, t) {
  return Dr(e) || Sr(e, t) || br(e, t) || yr();
}
function yr() {
  throw new TypeError(`Invalid attempt to destructure non-iterable instance.
In order to be iterable, non-array objects must have a [Symbol.iterator]() method.`);
}
function br(e, t) {
  if (e) {
    if (typeof e == "string")
      return We(e, t);
    var n = Object.prototype.toString.call(e).slice(8, -1);
    if (n === "Object" && e.constructor && (n = e.constructor.name), n === "Map" || n === "Set")
      return Array.from(e);
    if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n))
      return We(e, t);
  }
}
function We(e, t) {
  (t == null || t > e.length) && (t = e.length);
  for (var n = 0, r = new Array(t); n < t; n++)
    r[n] = e[n];
  return r;
}
function Sr(e, t) {
  var n = e == null ? null : typeof Symbol < "u" && e[Symbol.iterator] || e["@@iterator"];
  if (n != null) {
    var r = [], o = !0, i = !1, f, d;
    try {
      for (n = n.call(e); !(o = (f = n.next()).done) && (r.push(f.value), !(t && r.length === t)); o = !0)
        ;
    } catch (v) {
      i = !0, d = v;
    } finally {
      try {
        !o && n.return != null && n.return();
      } finally {
        if (i)
          throw d;
      }
    }
    return r;
  }
}
function Dr(e) {
  if (Array.isArray(e))
    return e;
}
function _r(e, t) {
  if (!(e instanceof t))
    throw new TypeError("Cannot call a class as a function");
}
function Ue(e, t) {
  for (var n = 0; n < t.length; n++) {
    var r = t[n];
    r.enumerable = r.enumerable || !1, r.configurable = !0, "value" in r && (r.writable = !0), Object.defineProperty(e, r.key, r);
  }
}
function wr(e, t, n) {
  return t && Ue(e.prototype, t), n && Ue(e, n), Object.defineProperty(e, "prototype", { writable: !1 }), e;
}
function Or(e, t) {
  if (typeof t != "function" && t !== null)
    throw new TypeError("Super expression must either be null or a function");
  e.prototype = Object.create(t && t.prototype, { constructor: { value: e, writable: !0, configurable: !0 } }), Object.defineProperty(e, "prototype", { writable: !1 }), t && Se(e, t);
}
function Se(e, t) {
  return Se = Object.setPrototypeOf || function(r, o) {
    return r.__proto__ = o, r;
  }, Se(e, t);
}
function Pr(e) {
  var t = xr();
  return function() {
    var r = se(e), o;
    if (t) {
      var i = se(this).constructor;
      o = Reflect.construct(r, arguments, i);
    } else
      o = r.apply(this, arguments);
    return Er(this, o);
  };
}
function Er(e, t) {
  if (t && (ue(t) === "object" || typeof t == "function"))
    return t;
  if (t !== void 0)
    throw new TypeError("Derived constructors may only return object or undefined");
  return C(e);
}
function C(e) {
  if (e === void 0)
    throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
  return e;
}
function xr() {
  if (typeof Reflect > "u" || !Reflect.construct || Reflect.construct.sham)
    return !1;
  if (typeof Proxy == "function")
    return !0;
  try {
    return Boolean.prototype.valueOf.call(Reflect.construct(Boolean, [], function() {
    })), !0;
  } catch {
    return !1;
  }
}
function se(e) {
  return se = Object.setPrototypeOf ? Object.getPrototypeOf : function(n) {
    return n.__proto__ || Object.getPrototypeOf(n);
  }, se(e);
}
function $(e, t, n) {
  return t in e ? Object.defineProperty(e, t, { value: n, enumerable: !0, configurable: !0, writable: !0 }) : e[t] = n, e;
}
var I = {
  touch: {
    start: "touchstart",
    move: "touchmove",
    stop: "touchend"
  },
  mouse: {
    start: "mousedown",
    move: "mousemove",
    stop: "mouseup"
  }
}, F = I.mouse, pe = /* @__PURE__ */ function(e) {
  Or(n, e);
  var t = Pr(n);
  function n() {
    var r;
    _r(this, n);
    for (var o = arguments.length, i = new Array(o), f = 0; f < o; f++)
      i[f] = arguments[f];
    return r = t.call.apply(t, [this].concat(i)), $(C(r), "state", {
      dragging: !1,
      // Used while dragging to determine deltas.
      lastX: NaN,
      lastY: NaN,
      touchIdentifier: null
    }), $(C(r), "mounted", !1), $(C(r), "handleDragStart", function(d) {
      if (r.props.onMouseDown(d), !r.props.allowAnyClick && typeof d.button == "number" && d.button !== 0)
        return !1;
      var v = r.findDOMNode();
      if (!v || !v.ownerDocument || !v.ownerDocument.body)
        throw new Error("<DraggableCore> not mounted on DragStart!");
      var b = v.ownerDocument;
      if (!(r.props.disabled || !(d.target instanceof b.defaultView.Node) || r.props.handle && !(0, x.matchesSelectorAndParentsTo)(d.target, r.props.handle, v) || r.props.cancel && (0, x.matchesSelectorAndParentsTo)(d.target, r.props.cancel, v))) {
        d.type === "touchstart" && d.preventDefault();
        var S = (0, x.getTouchIdentifier)(d);
        r.setState({
          touchIdentifier: S
        });
        var E = (0, X.getControlPosition)(d, S, C(r));
        if (E != null) {
          var O = E.x, j = E.y, A = (0, X.createCoreData)(C(r), O, j);
          (0, Z.default)("DraggableCore: handleDragStart: %j", A), (0, Z.default)("calling", r.props.onStart);
          var h = r.props.onStart(d, A);
          h === !1 || r.mounted === !1 || (r.props.enableUserSelectHack && (0, x.addUserSelectStyles)(b), r.setState({
            dragging: !0,
            lastX: O,
            lastY: j
          }), (0, x.addEvent)(b, F.move, r.handleDrag), (0, x.addEvent)(b, F.stop, r.handleDragStop));
        }
      }
    }), $(C(r), "handleDrag", function(d) {
      var v = (0, X.getControlPosition)(d, r.state.touchIdentifier, C(r));
      if (v != null) {
        var b = v.x, S = v.y;
        if (Array.isArray(r.props.grid)) {
          var E = b - r.state.lastX, O = S - r.state.lastY, j = (0, X.snapToGrid)(r.props.grid, E, O), A = Ie(j, 2);
          if (E = A[0], O = A[1], !E && !O)
            return;
          b = r.state.lastX + E, S = r.state.lastY + O;
        }
        var h = (0, X.createCoreData)(C(r), b, S);
        (0, Z.default)("DraggableCore: handleDrag: %j", h);
        var P = r.props.onDrag(d, h);
        if (P === !1 || r.mounted === !1) {
          try {
            r.handleDragStop(new MouseEvent("mouseup"));
          } catch {
            var D = document.createEvent("MouseEvents");
            D.initMouseEvent("mouseup", !0, !0, window, 0, 0, 0, 0, 0, !1, !1, !1, !1, 0, null), r.handleDragStop(D);
          }
          return;
        }
        r.setState({
          lastX: b,
          lastY: S
        });
      }
    }), $(C(r), "handleDragStop", function(d) {
      if (r.state.dragging) {
        var v = (0, X.getControlPosition)(d, r.state.touchIdentifier, C(r));
        if (v != null) {
          var b = v.x, S = v.y;
          if (Array.isArray(r.props.grid)) {
            var E = b - r.state.lastX || 0, O = S - r.state.lastY || 0, j = (0, X.snapToGrid)(r.props.grid, E, O), A = Ie(j, 2);
            E = A[0], O = A[1], b = r.state.lastX + E, S = r.state.lastY + O;
          }
          var h = (0, X.createCoreData)(C(r), b, S), P = r.props.onStop(d, h);
          if (P === !1 || r.mounted === !1)
            return !1;
          var D = r.findDOMNode();
          D && r.props.enableUserSelectHack && (0, x.removeUserSelectStyles)(D.ownerDocument), (0, Z.default)("DraggableCore: handleDragStop: %j", h), r.setState({
            dragging: !1,
            lastX: NaN,
            lastY: NaN
          }), D && ((0, Z.default)("DraggableCore: Removing handlers"), (0, x.removeEvent)(D.ownerDocument, F.move, r.handleDrag), (0, x.removeEvent)(D.ownerDocument, F.stop, r.handleDragStop));
        }
      }
    }), $(C(r), "onMouseDown", function(d) {
      return F = I.mouse, r.handleDragStart(d);
    }), $(C(r), "onMouseUp", function(d) {
      return F = I.mouse, r.handleDragStop(d);
    }), $(C(r), "onTouchStart", function(d) {
      return F = I.touch, r.handleDragStart(d);
    }), $(C(r), "onTouchEnd", function(d) {
      return F = I.touch, r.handleDragStop(d);
    }), r;
  }
  return wr(n, [{
    key: "componentDidMount",
    value: function() {
      this.mounted = !0;
      var o = this.findDOMNode();
      o && (0, x.addEvent)(o, I.touch.start, this.onTouchStart, {
        passive: !1
      });
    }
  }, {
    key: "componentWillUnmount",
    value: function() {
      this.mounted = !1;
      var o = this.findDOMNode();
      if (o) {
        var i = o.ownerDocument;
        (0, x.removeEvent)(i, I.mouse.move, this.handleDrag), (0, x.removeEvent)(i, I.touch.move, this.handleDrag), (0, x.removeEvent)(i, I.mouse.stop, this.handleDragStop), (0, x.removeEvent)(i, I.touch.stop, this.handleDragStop), (0, x.removeEvent)(o, I.touch.start, this.onTouchStart, {
          passive: !1
        }), this.props.enableUserSelectHack && (0, x.removeUserSelectStyles)(i);
      }
    }
    // React Strict Mode compatibility: if `nodeRef` is passed, we will use it instead of trying to find
    // the underlying DOM node ourselves. See the README for more information.
  }, {
    key: "findDOMNode",
    value: function() {
      var o, i, f;
      return (o = this.props) !== null && o !== void 0 && o.nodeRef ? (i = this.props) === null || i === void 0 || (f = i.nodeRef) === null || f === void 0 ? void 0 : f.current : hr.default.findDOMNode(this);
    }
  }, {
    key: "render",
    value: function() {
      return /* @__PURE__ */ ve.cloneElement(ve.Children.only(this.props.children), {
        // Note: mouseMove handler is attached to document so it will still function
        // when the user drags quickly and leaves the bounds of the element.
        onMouseDown: this.onMouseDown,
        onMouseUp: this.onMouseUp,
        // onTouchStart is added on `componentDidMount` so they can be added with
        // {passive: false}, which allows it to cancel. See
        // https://developers.google.com/web/updates/2017/01/scrolling-intervention
        onTouchEnd: this.onTouchEnd
      });
    }
  }]), n;
}(ve.Component);
fe.default = pe;
$(pe, "displayName", "DraggableCore");
$(pe, "propTypes", {
  /**
   * `allowAnyClick` allows dragging using any mouse button.
   * By default, we only accept the left button.
   *
   * Defaults to `false`.
   */
  allowAnyClick: T.default.bool,
  /**
   * `disabled`, if true, stops the <Draggable> from dragging. All handlers,
   * with the exception of `onMouseDown`, will not fire.
   */
  disabled: T.default.bool,
  /**
   * By default, we add 'user-select:none' attributes to the document body
   * to prevent ugly text selection during drag. If this is causing problems
   * for your app, set this to `false`.
   */
  enableUserSelectHack: T.default.bool,
  /**
   * `offsetParent`, if set, uses the passed DOM node to compute drag offsets
   * instead of using the parent node.
   */
  offsetParent: function(t, n) {
    if (t[n] && t[n].nodeType !== 1)
      throw new Error("Draggable's offsetParent must be a DOM Node.");
  },
  /**
   * `grid` specifies the x and y that dragging should snap to.
   */
  grid: T.default.arrayOf(T.default.number),
  /**
   * `handle` specifies a selector to be used as the handle that initiates drag.
   *
   * Example:
   *
   * ```jsx
   *   let App = React.createClass({
   *       render: function () {
   *         return (
   *            <Draggable handle=".handle">
   *              <div>
   *                  <div className="handle">Click me to drag</div>
   *                  <div>This is some other content</div>
   *              </div>
   *           </Draggable>
   *         );
   *       }
   *   });
   * ```
   */
  handle: T.default.string,
  /**
   * `cancel` specifies a selector to be used to prevent drag initialization.
   *
   * Example:
   *
   * ```jsx
   *   let App = React.createClass({
   *       render: function () {
   *           return(
   *               <Draggable cancel=".cancel">
   *                   <div>
   *                     <div className="cancel">You can't drag from here</div>
   *                     <div>Dragging here works fine</div>
   *                   </div>
   *               </Draggable>
   *           );
   *       }
   *   });
   * ```
   */
  cancel: T.default.string,
  /* If running in React Strict mode, ReactDOM.findDOMNode() is deprecated.
   * Unfortunately, in order for <Draggable> to work properly, we need raw access
   * to the underlying DOM node. If you want to avoid the warning, pass a `nodeRef`
   * as in this example:
   *
   * function MyComponent() {
   *   const nodeRef = React.useRef(null);
   *   return (
   *     <Draggable nodeRef={nodeRef}>
   *       <div ref={nodeRef}>Example Target</div>
   *     </Draggable>
   *   );
   * }
   *
   * This can be used for arbitrarily nested components, so long as the ref ends up
   * pointing to the actual child DOM node and not a custom component.
   */
  nodeRef: T.default.object,
  /**
   * Called when dragging starts.
   * If this function returns the boolean false, dragging will be canceled.
   */
  onStart: T.default.func,
  /**
   * Called while dragging.
   * If this function returns the boolean false, dragging will be canceled.
   */
  onDrag: T.default.func,
  /**
   * Called when dragging stops.
   * If this function returns the boolean false, the drag will remain active.
   */
  onStop: T.default.func,
  /**
   * A workaround option which can be passed if onMouseDown needs to be accessed,
   * since it'll always be blocked (as there is internal use of onMouseDown)
   */
  onMouseDown: T.default.func,
  /**
   * `scale`, if set, applies scaling while dragging an element
   */
  scale: T.default.number,
  /**
   * These properties should be defined on the child, not here.
   */
  className: ye.dontSetMe,
  style: ye.dontSetMe,
  transform: ye.dontSetMe
});
$(pe, "defaultProps", {
  allowAnyClick: !1,
  // by default only accept left click
  disabled: !1,
  enableUserSelectHack: !0,
  onStart: function() {
  },
  onDrag: function() {
  },
  onStop: function() {
  },
  onMouseDown: function() {
  },
  scale: 1
});
(function(e) {
  function t(a) {
    "@babel/helpers - typeof";
    return t = typeof Symbol == "function" && typeof Symbol.iterator == "symbol" ? function(l) {
      return typeof l;
    } : function(l) {
      return l && typeof Symbol == "function" && l.constructor === Symbol && l !== Symbol.prototype ? "symbol" : typeof l;
    }, t(a);
  }
  Object.defineProperty(e, "__esModule", {
    value: !0
  }), Object.defineProperty(e, "DraggableCore", {
    enumerable: !0,
    get: function() {
      return b.default;
    }
  }), e.default = void 0;
  var n = A(p), r = O(Fe), o = O(qe), i = O(It), f = _, d = U, v = W, b = O(fe), S = O(de), E = ["axis", "bounds", "children", "defaultPosition", "defaultClassName", "defaultClassNameDragging", "defaultClassNameDragged", "position", "positionOffset", "scale"];
  function O(a) {
    return a && a.__esModule ? a : { default: a };
  }
  function j(a) {
    if (typeof WeakMap != "function")
      return null;
    var l = /* @__PURE__ */ new WeakMap(), u = /* @__PURE__ */ new WeakMap();
    return (j = function(c) {
      return c ? u : l;
    })(a);
  }
  function A(a, l) {
    if (!l && a && a.__esModule)
      return a;
    if (a === null || t(a) !== "object" && typeof a != "function")
      return { default: a };
    var u = j(l);
    if (u && u.has(a))
      return u.get(a);
    var s = {}, c = Object.defineProperty && Object.getOwnPropertyDescriptor;
    for (var g in a)
      if (g !== "default" && Object.prototype.hasOwnProperty.call(a, g)) {
        var m = c ? Object.getOwnPropertyDescriptor(a, g) : null;
        m && (m.get || m.set) ? Object.defineProperty(s, g, m) : s[g] = a[g];
      }
    return s.default = a, u && u.set(a, s), s;
  }
  function h() {
    return h = Object.assign || function(a) {
      for (var l = 1; l < arguments.length; l++) {
        var u = arguments[l];
        for (var s in u)
          Object.prototype.hasOwnProperty.call(u, s) && (a[s] = u[s]);
      }
      return a;
    }, h.apply(this, arguments);
  }
  function P(a, l) {
    if (a == null)
      return {};
    var u = D(a, l), s, c;
    if (Object.getOwnPropertySymbols) {
      var g = Object.getOwnPropertySymbols(a);
      for (c = 0; c < g.length; c++)
        s = g[c], !(l.indexOf(s) >= 0) && Object.prototype.propertyIsEnumerable.call(a, s) && (u[s] = a[s]);
    }
    return u;
  }
  function D(a, l) {
    if (a == null)
      return {};
    var u = {}, s = Object.keys(a), c, g;
    for (g = 0; g < s.length; g++)
      c = s[g], !(l.indexOf(c) >= 0) && (u[c] = a[c]);
    return u;
  }
  function k(a, l) {
    var u = Object.keys(a);
    if (Object.getOwnPropertySymbols) {
      var s = Object.getOwnPropertySymbols(a);
      l && (s = s.filter(function(c) {
        return Object.getOwnPropertyDescriptor(a, c).enumerable;
      })), u.push.apply(u, s);
    }
    return u;
  }
  function M(a) {
    for (var l = 1; l < arguments.length; l++) {
      var u = arguments[l] != null ? arguments[l] : {};
      l % 2 ? k(Object(u), !0).forEach(function(s) {
        L(a, s, u[s]);
      }) : Object.getOwnPropertyDescriptors ? Object.defineProperties(a, Object.getOwnPropertyDescriptors(u)) : k(Object(u)).forEach(function(s) {
        Object.defineProperty(a, s, Object.getOwnPropertyDescriptor(u, s));
      });
    }
    return a;
  }
  function J(a, l) {
    return st(a) || ut(a, l) || re(a, l) || te();
  }
  function te() {
    throw new TypeError(`Invalid attempt to destructure non-iterable instance.
In order to be iterable, non-array objects must have a [Symbol.iterator]() method.`);
  }
  function re(a, l) {
    if (a) {
      if (typeof a == "string")
        return Te(a, l);
      var u = Object.prototype.toString.call(a).slice(8, -1);
      if (u === "Object" && a.constructor && (u = a.constructor.name), u === "Map" || u === "Set")
        return Array.from(a);
      if (u === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(u))
        return Te(a, l);
    }
  }
  function Te(a, l) {
    (l == null || l > a.length) && (l = a.length);
    for (var u = 0, s = new Array(l); u < l; u++)
      s[u] = a[u];
    return s;
  }
  function ut(a, l) {
    var u = a == null ? null : typeof Symbol < "u" && a[Symbol.iterator] || a["@@iterator"];
    if (u != null) {
      var s = [], c = !0, g = !1, m, w;
      try {
        for (u = u.call(a); !(c = (m = u.next()).done) && (s.push(m.value), !(l && s.length === l)); c = !0)
          ;
      } catch (y) {
        g = !0, w = y;
      } finally {
        try {
          !c && u.return != null && u.return();
        } finally {
          if (g)
            throw w;
        }
      }
      return s;
    }
  }
  function st(a) {
    if (Array.isArray(a))
      return a;
  }
  function ct(a, l) {
    if (!(a instanceof l))
      throw new TypeError("Cannot call a class as a function");
  }
  function Ne(a, l) {
    for (var u = 0; u < l.length; u++) {
      var s = l[u];
      s.enumerable = s.enumerable || !1, s.configurable = !0, "value" in s && (s.writable = !0), Object.defineProperty(a, s.key, s);
    }
  }
  function ft(a, l, u) {
    return l && Ne(a.prototype, l), u && Ne(a, u), Object.defineProperty(a, "prototype", { writable: !1 }), a;
  }
  function dt(a, l) {
    if (typeof l != "function" && l !== null)
      throw new TypeError("Super expression must either be null or a function");
    a.prototype = Object.create(l && l.prototype, { constructor: { value: a, writable: !0, configurable: !0 } }), Object.defineProperty(a, "prototype", { writable: !1 }), l && ge(a, l);
  }
  function ge(a, l) {
    return ge = Object.setPrototypeOf || function(s, c) {
      return s.__proto__ = c, s;
    }, ge(a, l);
  }
  function pt(a) {
    var l = mt();
    return function() {
      var s = ne(a), c;
      if (l) {
        var g = ne(this).constructor;
        c = Reflect.construct(s, arguments, g);
      } else
        c = s.apply(this, arguments);
      return gt(this, c);
    };
  }
  function gt(a, l) {
    if (l && (t(l) === "object" || typeof l == "function"))
      return l;
    if (l !== void 0)
      throw new TypeError("Derived constructors may only return object or undefined");
    return B(a);
  }
  function B(a) {
    if (a === void 0)
      throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
    return a;
  }
  function mt() {
    if (typeof Reflect > "u" || !Reflect.construct || Reflect.construct.sham)
      return !1;
    if (typeof Proxy == "function")
      return !0;
    try {
      return Boolean.prototype.valueOf.call(Reflect.construct(Boolean, [], function() {
      })), !0;
    } catch {
      return !1;
    }
  }
  function ne(a) {
    return ne = Object.setPrototypeOf ? Object.getPrototypeOf : function(u) {
      return u.__proto__ || Object.getPrototypeOf(u);
    }, ne(a);
  }
  function L(a, l, u) {
    return l in a ? Object.defineProperty(a, l, { value: u, enumerable: !0, configurable: !0, writable: !0 }) : a[l] = u, a;
  }
  var ae = /* @__PURE__ */ function(a) {
    dt(u, a);
    var l = pt(u);
    function u(s) {
      var c;
      return ct(this, u), c = l.call(this, s), L(B(c), "onDragStart", function(g, m) {
        (0, S.default)("Draggable: onDragStart: %j", m);
        var w = c.props.onStart(g, (0, d.createDraggableData)(B(c), m));
        if (w === !1)
          return !1;
        c.setState({
          dragging: !0,
          dragged: !0
        });
      }), L(B(c), "onDrag", function(g, m) {
        if (!c.state.dragging)
          return !1;
        (0, S.default)("Draggable: onDrag: %j", m);
        var w = (0, d.createDraggableData)(B(c), m), y = {
          x: w.x,
          y: w.y
        };
        if (c.props.bounds) {
          var Q = y.x, G = y.y;
          y.x += c.state.slackX, y.y += c.state.slackY;
          var K = (0, d.getBoundPosition)(B(c), y.x, y.y), H = J(K, 2), me = H[0], oe = H[1];
          y.x = me, y.y = oe, y.slackX = c.state.slackX + (Q - y.x), y.slackY = c.state.slackY + (G - y.y), w.x = y.x, w.y = y.y, w.deltaX = y.x - c.state.x, w.deltaY = y.y - c.state.y;
        }
        var ie = c.props.onDrag(g, w);
        if (ie === !1)
          return !1;
        c.setState(y);
      }), L(B(c), "onDragStop", function(g, m) {
        if (!c.state.dragging)
          return !1;
        var w = c.props.onStop(g, (0, d.createDraggableData)(B(c), m));
        if (w === !1)
          return !1;
        (0, S.default)("Draggable: onDragStop: %j", m);
        var y = {
          dragging: !1,
          slackX: 0,
          slackY: 0
        }, Q = !!c.props.position;
        if (Q) {
          var G = c.props.position, K = G.x, H = G.y;
          y.x = K, y.y = H;
        }
        c.setState(y);
      }), c.state = {
        // Whether or not we are currently dragging.
        dragging: !1,
        // Whether or not we have been dragged before.
        dragged: !1,
        // Current transform x and y.
        x: s.position ? s.position.x : s.defaultPosition.x,
        y: s.position ? s.position.y : s.defaultPosition.y,
        prevPropsPosition: M({}, s.position),
        // Used for compensating for out-of-bounds drags
        slackX: 0,
        slackY: 0,
        // Can only determine if SVG after mounting
        isElementSVG: !1
      }, s.position && !(s.onDrag || s.onStop) && console.warn("A `position` was applied to this <Draggable>, without drag handlers. This will make this component effectively undraggable. Please attach `onDrag` or `onStop` handlers so you can adjust the `position` of this element."), c;
    }
    return ft(u, [{
      key: "componentDidMount",
      value: function() {
        typeof window.SVGElement < "u" && this.findDOMNode() instanceof window.SVGElement && this.setState({
          isElementSVG: !0
        });
      }
    }, {
      key: "componentWillUnmount",
      value: function() {
        this.setState({
          dragging: !1
        });
      }
      // React Strict Mode compatibility: if `nodeRef` is passed, we will use it instead of trying to find
      // the underlying DOM node ourselves. See the README for more information.
    }, {
      key: "findDOMNode",
      value: function() {
        var c, g, m;
        return (c = (g = this.props) === null || g === void 0 || (m = g.nodeRef) === null || m === void 0 ? void 0 : m.current) !== null && c !== void 0 ? c : o.default.findDOMNode(this);
      }
    }, {
      key: "render",
      value: function() {
        var c, g = this.props;
        g.axis, g.bounds;
        var m = g.children, w = g.defaultPosition, y = g.defaultClassName, Q = g.defaultClassNameDragging, G = g.defaultClassNameDragged, K = g.position, H = g.positionOffset;
        g.scale;
        var me = P(g, E), oe = {}, ie = null, ht = !!K, Re = !ht || this.state.dragging, je = K || w, Ae = {
          // Set left if horizontal drag is enabled
          x: (0, d.canDragX)(this) && Re ? this.state.x : je.x,
          // Set top if vertical drag is enabled
          y: (0, d.canDragY)(this) && Re ? this.state.y : je.y
        };
        this.state.isElementSVG ? ie = (0, f.createSVGTransform)(Ae, H) : oe = (0, f.createCSSTransform)(Ae, H);
        var vt = (0, i.default)(m.props.className || "", y, (c = {}, L(c, Q, this.state.dragging), L(c, G, this.state.dragged), c));
        return /* @__PURE__ */ n.createElement(b.default, h({}, me, {
          onStart: this.onDragStart,
          onDrag: this.onDrag,
          onStop: this.onDragStop
        }), /* @__PURE__ */ n.cloneElement(n.Children.only(m), {
          className: vt,
          style: M(M({}, m.props.style), oe),
          transform: ie
        }));
      }
    }], [{
      key: "getDerivedStateFromProps",
      value: (
        // React 16.3+
        // Arity (props, state)
        function(c, g) {
          var m = c.position, w = g.prevPropsPosition;
          return m && (!w || m.x !== w.x || m.y !== w.y) ? ((0, S.default)("Draggable: getDerivedStateFromProps %j", {
            position: m,
            prevPropsPosition: w
          }), {
            x: m.x,
            y: m.y,
            prevPropsPosition: M({}, m)
          }) : null;
        }
      )
    }]), u;
  }(n.Component);
  e.default = ae, L(ae, "displayName", "Draggable"), L(ae, "propTypes", M(M({}, b.default.propTypes), {}, {
    /**
     * `axis` determines which axis the draggable can move.
     *
     *  Note that all callbacks will still return data as normal. This only
     *  controls flushing to the DOM.
     *
     * 'both' allows movement horizontally and vertically.
     * 'x' limits movement to horizontal axis.
     * 'y' limits movement to vertical axis.
     * 'none' limits all movement.
     *
     * Defaults to 'both'.
     */
    axis: r.default.oneOf(["both", "x", "y", "none"]),
    /**
     * `bounds` determines the range of movement available to the element.
     * Available values are:
     *
     * 'parent' restricts movement within the Draggable's parent node.
     *
     * Alternatively, pass an object with the following properties, all of which are optional:
     *
     * {left: LEFT_BOUND, right: RIGHT_BOUND, bottom: BOTTOM_BOUND, top: TOP_BOUND}
     *
     * All values are in px.
     *
     * Example:
     *
     * ```jsx
     *   let App = React.createClass({
     *       render: function () {
     *         return (
     *            <Draggable bounds={{right: 300, bottom: 300}}>
     *              <div>Content</div>
     *           </Draggable>
     *         );
     *       }
     *   });
     * ```
     */
    bounds: r.default.oneOfType([r.default.shape({
      left: r.default.number,
      right: r.default.number,
      top: r.default.number,
      bottom: r.default.number
    }), r.default.string, r.default.oneOf([!1])]),
    defaultClassName: r.default.string,
    defaultClassNameDragging: r.default.string,
    defaultClassNameDragged: r.default.string,
    /**
     * `defaultPosition` specifies the x and y that the dragged item should start at
     *
     * Example:
     *
     * ```jsx
     *      let App = React.createClass({
     *          render: function () {
     *              return (
     *                  <Draggable defaultPosition={{x: 25, y: 25}}>
     *                      <div>I start with transformX: 25px and transformY: 25px;</div>
     *                  </Draggable>
     *              );
     *          }
     *      });
     * ```
     */
    defaultPosition: r.default.shape({
      x: r.default.number,
      y: r.default.number
    }),
    positionOffset: r.default.shape({
      x: r.default.oneOfType([r.default.number, r.default.string]),
      y: r.default.oneOfType([r.default.number, r.default.string])
    }),
    /**
     * `position`, if present, defines the current position of the element.
     *
     *  This is similar to how form elements in React work - if no `position` is supplied, the component
     *  is uncontrolled.
     *
     * Example:
     *
     * ```jsx
     *      let App = React.createClass({
     *          render: function () {
     *              return (
     *                  <Draggable position={{x: 25, y: 25}}>
     *                      <div>I start with transformX: 25px and transformY: 25px;</div>
     *                  </Draggable>
     *              );
     *          }
     *      });
     * ```
     */
    position: r.default.shape({
      x: r.default.number,
      y: r.default.number
    }),
    /**
     * These properties should be defined on the child, not here.
     */
    className: v.dontSetMe,
    style: v.dontSetMe,
    transform: v.dontSetMe
  })), L(ae, "defaultProps", M(M({}, b.default.defaultProps), {}, {
    axis: "both",
    bounds: !1,
    defaultClassName: "react-draggable",
    defaultClassNameDragging: "react-draggable-dragging",
    defaultClassNameDragged: "react-draggable-dragged",
    defaultPosition: {
      x: 0,
      y: 0
    },
    scale: 1
  }));
})(Ve);
var nt = Ve, at = nt.default, Cr = nt.DraggableCore;
ce.exports = at;
ce.exports.default = at;
ce.exports.DraggableCore = Cr;
var Mr = ce.exports;
const Tr = /* @__PURE__ */ St(Mr), Nr = He()((e) => ({
  closeButton: {
    position: "absolute",
    right: e.spacing(1),
    top: e.spacing(1),
    color: e.palette.grey[500]
  }
}));
function Rr(e) {
  return p.createElement(
    Tr,
    { handle: "#draggable-dialog-title", cancel: '[class*="MuiDialogContent-root"]' },
    p.createElement(Ct, { ...e })
  );
}
function jr(e) {
  const { classes: t } = Nr(), { title: n, children: r, onClose: o } = e;
  return p.createElement(
    _t,
    {
      ...e,
      PaperComponent: Rr,
      "aria-labelledby": "draggable-dialog-title"
      // this area is important for the draggable functionality
    },
    p.createElement(
      wt,
      null,
      p.createElement(
        Ot,
        { style: { cursor: "move" }, id: "draggable-dialog-title" },
        n,
        o ? p.createElement(
          Pt,
          { className: t.closeButton, onClick: () => {
            o();
          } },
          p.createElement(Et, null)
        ) : null
      ),
      p.createElement(xt, null),
      r
    )
  );
}
const Ar = Dt(jr);
function Le(e, t, n = 1) {
  const r = t.map((i) => e.findIndex((f) => f.name === i)).sort((i, f) => i - f);
  let o = 0;
  for (let i = 0; i < r.length; i++) {
    const f = r[i], d = Math.max(o, f - n);
    d >= o && e.splice(d, 0, e.splice(f, 1)[0]), o = o + 1;
  }
  return e;
}
function Ye(e, t, n = 1) {
  const r = t.map((i) => e.findIndex((f) => f.name === i)).sort((i, f) => f - i);
  let o = e.length - 1;
  for (let i = 0; i < r.length; i++) {
    const f = r[i], d = Math.min(o, f + n);
    d <= o && e.splice(d, 0, e.splice(f, 1)[0]), o = o - 1;
  }
  return e;
}
var xe = {}, kr = _e;
Object.defineProperty(xe, "__esModule", {
  value: !0
});
var ot = xe.default = void 0, $r = kr(De()), Be = we, Ir = (0, $r.default)([/* @__PURE__ */ (0, Be.jsx)("path", {
  d: "M6 17.59 7.41 19 12 14.42 16.59 19 18 17.59l-6-6z"
}, "0"), /* @__PURE__ */ (0, Be.jsx)("path", {
  d: "m6 11 1.41 1.41L12 7.83l4.59 4.58L18 11l-6-6z"
}, "1")], "KeyboardDoubleArrowUp");
ot = xe.default = Ir;
var Ce = {}, Wr = _e;
Object.defineProperty(Ce, "__esModule", {
  value: !0
});
var it = Ce.default = void 0, Ur = Wr(De()), Xe = we, Lr = (0, Ur.default)([/* @__PURE__ */ (0, Xe.jsx)("path", {
  d: "M18 6.41 16.59 5 12 9.58 7.41 5 6 6.41l6 6z"
}, "0"), /* @__PURE__ */ (0, Xe.jsx)("path", {
  d: "m18 13-1.41-1.41L12 16.17l-4.59-4.58L6 13l6 6z"
}, "1")], "KeyboardDoubleArrowDown");
it = Ce.default = Lr;
var Me = {}, Yr = _e;
Object.defineProperty(Me, "__esModule", {
  value: !0
});
var lt = Me.default = void 0, Br = Yr(De()), Xr = we, Fr = (0, Br.default)(/* @__PURE__ */ (0, Xr.jsx)("path", {
  d: "M7.41 15.41 12 10.83l4.59 4.58L18 14l-6-6-6 6z"
}), "KeyboardArrowUp");
lt = Me.default = Fr;
const Hr = He()({
  content: {
    minWidth: 800
  }
});
function Qr({ model: e, handleClose: t }) {
  const { classes: n } = Hr(), { sources: r } = e, [o, i] = ee(Mt(r || [])), [f, d] = q.useLocalStorage("multiwiggle-showTips", !0);
  return p.createElement(
    Ar,
    { open: !0, onClose: t, maxWidth: "xl", title: "Multi-wiggle color/arrangement editor" },
    p.createElement(
      Tt,
      { className: n.content },
      p.createElement(Y, { variant: "contained", style: { float: "right" }, onClick: () => d(!f) }, f ? "Hide tips" : "Show tips"),
      p.createElement("br", null),
      f ? p.createElement(
        p.Fragment,
        null,
        "Helpful tips",
        p.createElement(
          "ul",
          null,
          p.createElement("li", null, "You can select rows in the table with the checkboxes"),
          p.createElement("li", null, "Multi-select is enabled with shift-click and control-click"),
          p.createElement("li", null, 'The "Move selected items up/down" can re-arrange subtracks'),
          p.createElement("li", null, "Sorting the data grid itself can also re-arrange subtracks"),
          p.createElement("li", null, "Changes are applied when you hit Submit"),
          p.createElement("li", null, "You can click and drag the dialog box to move it on the screen"),
          p.createElement("li", null, "Columns in the table can be hidden using a vertical '...' menu on the right side of each column")
        )
      ) : null,
      p.createElement(qr, { rows: o, onChange: i, showTips: f })
    ),
    p.createElement(
      Nt,
      null,
      p.createElement(Y, { variant: "contained", type: "submit", color: "inherit", onClick: () => {
        e.clearLayout(), i(e.sources || []);
      } }, "Clear custom settings"),
      p.createElement(Y, { variant: "contained", color: "secondary", onClick: () => {
        t(), i([...e.sources || []]);
      } }, "Cancel"),
      p.createElement(Y, { variant: "contained", color: "primary", type: "submit", onClick: () => {
        e.setLayout(o), t();
      } }, "Submit")
    )
  );
}
function qr({ rows: e, onChange: t, showTips: n }) {
  const [r, o] = ee(null), [i, f] = ee([]), { name: d, color: v, baseUri: b, ...S } = e[0], [E, O] = ee("blue"), [j, A] = ee({
    idx: 0,
    field: null
  });
  return p.createElement(
    "div",
    null,
    p.createElement(Y, { disabled: !i.length, onClick: (h) => o(h.currentTarget) }, "Change color of selected items"),
    p.createElement(
      Y,
      { onClick: () => t(Le([...e], i)), disabled: !i.length },
      p.createElement(lt, null),
      n ? "Move selected items up" : null
    ),
    p.createElement(
      Y,
      { onClick: () => t(Ye([...e], i)), disabled: !i.length },
      p.createElement(Rt, null),
      n ? "Move selected items down" : null
    ),
    p.createElement(
      Y,
      { onClick: () => t(Le([...e], i, e.length)), disabled: !i.length },
      p.createElement(ot, null),
      n ? "Move selected items to top" : null
    ),
    p.createElement(
      Y,
      { onClick: () => t(Ye([...e], i, e.length)), disabled: !i.length },
      p.createElement(it, null),
      n ? "Move selected items to bottom" : null
    ),
    p.createElement(kt.ColorPopover, { anchorEl: r, color: E, onChange: (h) => {
      O(h), i.forEach((P) => {
        const D = e.find((k) => k.name === P);
        D && (D.color = h);
      }), t([...e]);
    }, onClose: () => o(null) }),
    p.createElement(
      "div",
      { style: { height: 400, width: "100%" } },
      p.createElement(jt, { getRowId: (h) => h.name, checkboxSelection: !0, disableRowSelectionOnClick: !0, onRowSelectionModelChange: (h) => f(h), rows: e, rowHeight: 25, columnHeaderHeight: 33, columns: [
        {
          field: "color",
          headerName: "Color",
          renderCell: (h) => {
            const { value: P, id: D } = h;
            return p.createElement($t, { color: P || "blue", onChange: (k) => {
              const M = e.find((J) => J.name === D);
              M && (M.color = k), t([...e]);
            } });
          }
        },
        {
          field: "name",
          sortingOrder: [null],
          headerName: "Name",
          width: q.measureGridWidth(e.map((h) => h.name))
        },
        ...Object.keys(S).map((h) => ({
          field: h,
          sortingOrder: [null],
          renderCell: (P) => {
            const { value: D } = P;
            return q.isUriLocation(D) ? p.createElement(At.UriLink, { value: D }) : q.getStr(D);
          },
          // @ts-ignore
          width: q.measureGridWidth(e.map((P) => P[h]))
        }))
      ], sortModel: [
        /* we control the sort as a controlled component using onSortModelChange */
      ], onSortModelChange: (h) => {
        const P = h[0], D = (j.idx + 1) % 2, k = (P == null ? void 0 : P.field) || j.field;
        A({ idx: D, field: k }), t(k ? [...e].sort((M, J) => {
          const te = q.getStr(M[k]), re = q.getStr(J[k]);
          return D === 1 ? te.localeCompare(re) : re.localeCompare(te);
        }) : e);
      } })
    )
  );
}
export {
  Qr as default
};
