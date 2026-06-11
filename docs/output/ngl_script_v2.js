document.addEventListener('DOMContentLoaded', function() {

  console.log('✅ NGL script running');

  var stage = new NGL.Stage('ngl_viewer', {
    backgroundColor: 'white'
  });

  var component = null;
  window.current_receptor = 'FcgR1_fucos';

  // -----------------------------
  // LOAD STRUCTURE ONCE ONLY ✅
  // -----------------------------
  stage.loadFile('data/Fc_structure.pdb', { ext: 'pdb' })
    .then(function(comp) {

      console.log('✅ Structure loaded');

      component = comp;

      // Cartoon backbone
      comp.addRepresentation('cartoon', {
        colorScheme: 'chainname'
      });

      // Base surface (very light)
      comp.addRepresentation('surface', {
        sele: ':A or :B',
        opacity: 0.01
      });

      comp.autoView();

      // ✅ Initial highlight
      applyHighlights('FcgR1_fucos');
    })
    .catch(function(err) {
      console.error('❌ Failed to load structure:', err);
    });


  // -----------------------------
// ✅ store highlight representations globally
var highlightReps = [];

// -----------------------------
// HIGHLIGHT FUNCTION ✅
// -----------------------------
function applyHighlights(receptor) {

  console.log('Available highlight keys:', Object.keys(window.highlight_sites));
  console.log('Received receptor:', receptor);

  if (!component) return;

  // ✅ FIXED key mapping (robust)
  var keys = Object.keys(window.highlight_sites);

  var key;
  if (keys.includes(receptor)) {
    key = receptor;
  } else {
    key = keys.find(k => k !== 'FcgR1_fucos');
  }

  console.log('Using key:', key);

  var sites = window.highlight_sites[key];
  if (!sites) return;

// ✅ ensure numeric residue IDs
var blueList = (sites.blue || []).map(Number).filter(n => !isNaN(n));
var redList  = (sites.red  || []).map(Number).filter(n => !isNaN(n));

// ✅ build proper NGL selection strings
var blue = blueList.map(n => 'resi ' + n).join(' or ');
var red  = redList.map(n => 'resi ' + n).join(' or ');

console.log('🔵 blueList:', blueList);
console.log('🔴 redList:', redList);
console.log('🔵 blue selection:', blue);
console.log('🔴 red selection:', red);

  // ✅ PROPERLY remove previous highlight reps
  highlightReps.forEach(function(rep) {
    component.removeRepresentation(rep);
  });
  highlightReps = [];

  // ✅ add blue highlights
  if (blue.length > 0) {
    var rep1 = component.addRepresentation('surface', {
      sele: '(' + blue + ') and (:A or :B)',
      color: '#4F8EDB',
      opacity: 0.15
    });
    highlightReps.push(rep1);
  }

  // ✅ add red highlights
  if (red.length > 0) {
    var rep2 = component.addRepresentation('surface', {
      sele: '(' + red + ') and (:A or :B)',
      color: '#E36A6A',
      opacity: 0.15
    });
    highlightReps.push(rep2);
  }
}


  // -----------------------------
  // DROPDOWN → STRUCTURE ✅
  // -----------------------------
  var plot = document.querySelector('.js-plotly-plot');

  if (plot) {

    plot.on('plotly_buttonclicked', function(ev) {

      if (!ev || !ev.button || !ev.button.label) return;

      var label = ev.button.label;

      console.log('🔄 Dropdown:', label);

      window.current_receptor = label;
      applyHighlights(label);
    });
  }


  // -----------------------------
  // CLICK → ZOOM + HIGHLIGHT ✅
  // -----------------------------
  if (plot) {

    plot.on('plotly_click', function(data) {

      if (!data.points || data.points.length === 0) return;

      var site = data.points[0].x;
      console.log('🟡 Clicked site:', site);

      var sele = site + ' and (:A or :B)';

      component.addRepresentation('spacefill', {
        sele: sele,
        color: 'green',
        radius: 2.0
      });

      component.autoView(sele);
    });
  }

});
