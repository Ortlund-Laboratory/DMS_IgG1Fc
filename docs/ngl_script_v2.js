document.addEventListener('DOMContentLoaded', function() {

  console.log('✅ NGL script running');

  var stage = new NGL.Stage('ngl_viewer', {
    backgroundColor: 'white'
  });

  window.component = null;
  window.current_receptor = 'FcgR1_fucos';

  // ✅ store highlight representations globally
  window.highlightReps = [];

  // -----------------------------
  // LOAD STRUCTURE ✅
  // -----------------------------
  stage.loadFile('data/Fc_structure.pdb', { ext: 'pdb' })
    .then(function(comp) {

      console.log('✅ Structure loaded');

      window.component = comp;

      // Cartoon backbone
// Base (keep yellowish wheat tone)
comp.addRepresentation('cartoon', {
  colorScheme: 'uniform',
  color: '#c4a484'
});

// Dark wheat
comp.addRepresentation('cartoon', {
  sele: ':A',
  color: '#c2a76d'
});

// Light wheat
comp.addRepresentation('cartoon', {
  sele: ':B',
  color: '#f5deb3'
});

      // Base surface
      comp.addRepresentation('surface', {
        sele: ':A or :B',
        opacity: 0.01
      });

      comp.autoView();



    })
    .catch(function(err) {
      console.error('❌ Failed to load structure:', err);
    });




  // -----------------------------
  // HIGHLIGHT FUNCTION ✅
  // -----------------------------
/*  
function applyHighlights(receptor) {

    console.log('Available highlight keys:', Object.keys(window.highlight_sites));
    console.log('Received receptor:', receptor);

    if (!window.component) return;

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

    var blueList = (sites.blue || []).map(Number).filter(n => !isNaN(n));
    var redList  = (sites.red  || []).map(Number).filter(n => !isNaN(n));

    var blue = blueList.map(n => 'resi ' + n).join(' or ');
    var red  = redList.map(n => 'resi ' + n).join(' or ');

    console.log('🔵 blueList:', blueList);
    console.log('🔴 redList:', redList);

    // ✅ remove previous highlights
    window.highlightReps = [];

    // ✅ blue highlights
    if (blue.length > 0) {
      var rep1 = window.component.addRepresentation('surface', {
        sele: '(' + blue + ') and (:A or :B)',
        color: '#4F8EDB',
        opacity: 0.15
      });
      window.highlightReps.push(rep1);
    }

    // ✅ red highlights
    if (red.length > 0) {
      var rep2 = window.component.addRepresentation('surface', {
        sele: '(' + red + ') and (:A or :B)',
        color: '#E36A6A',
        opacity: 0.15
      });
      window.highlightReps.push(rep2);
    }
  }
*/
});

window.highlightFromPlot = function(site) {

  if (!window.component) {
    console.log("❌ component not ready");
    return;
  }

  var resid = parseInt(site);

  var sele = 'resi ' + resid;

  console.log("🎯 NGL selecting:", sele);

  // ✅ remove previous click highlight
  if (window.lastClickRep) {
    window.component.removeRepresentation(window.lastClickRep);
  }

  // ✅ add new highlight (green)
  window.lastClickRep = window.component.addRepresentation('ball+stick', {
    sele: sele,
    color: 'green',
    scale: 4.0,
    aspectRatio: 2.0,
    radius: 0.2
  });

  // ✅ zoom to this residue ONLY
  window.component.autoView(sele);

};

window.updateStructureHighlights = function({ blue, red }) {

  console.log("Updating structure", blue.length, red.length);

  if (!window.component) {
    console.warn("Structure not ready");
    return;
  }

  const seleBlue = blue.map(n => 'resi ' + n).join(' or ');
  const seleRed  = red.map(n => 'resi ' + n).join(' or ');

  // ✅ Remove only previous highlight layers
  if (window.highlightReps && window.highlightReps.length > 0) {
    window.highlightReps.forEach(rep => window.component.removeRepresentation(rep));
  }

  window.highlightReps = [];

  // ✅ Blue highlights
  if (blue.length > 0) {
    const rep1 = window.component.addRepresentation('surface', {
      sele: '(' + seleBlue + ') and (:A or :B)',
      color: '#4F8EDB',
      opacity: 0.15
    });
    window.highlightReps.push(rep1);
  }

  // ✅ Red highlights
  if (red.length > 0) {
    const rep2 = window.component.addRepresentation('surface', {
      sele: '(' + seleRed + ') and (:A or :B)',
      color: '#E36A6A',
      opacity: 0.15
    });
    window.highlightReps.push(rep2);
  };
};
