document.addEventListener('DOMContentLoaded', function() {

  console.log('✅ NGL script running:', typeof NGL);

  var stage = new NGL.Stage('ngl_viewer', { backgroundColor: 'white' });

  function loadStructure(receptor) {

    console.log('🔁 Loading structure for:', receptor);

    stage.loadFile('data/Fc_structure.pdb', { ext: 'pdb' })
      .then(function(comp) {

        console.log('✅ Structure loaded into NGL');

        comp.addRepresentation('cartoon', {
          colorScheme: 'chainname'
        });

        comp.addRepresentation('surface', {
          sele: ':A or :B',
          opacity: 0.01
        });

        // ✅ HARD FILTER: only allow valid receptors
        if (!(receptor in window.highlight_sites)) {
         console.log('⛔ Ignoring invalid receptor call:', receptor);
         return;
       }

       var sites = window.highlight_sites[receptor];

	var blue = sites.blue.join(' ');
        var red  = sites.red.join(' ');
        
	if (blue.length > 0) {
          comp.addRepresentation('surface', {
            sele: '(' + blue + ') and (:A or :B)',
            color: '#4F8EDB',
            opacity: 0.05
          });
        }

        if (red.length > 0) {
          comp.addRepresentation('surface', {
            sele: '(' + red + ') and (:A or :B)',
            color: '#E36A6A',
            opacity: 0.05
          });
        }

        comp.autoView();
      })
      .catch(function(err) {
        console.log('❌ Failed to load structure:', err);
      });
  }

  loadStructure('FcgR2b');

	// ✅ CONNECT HEATMAP CLICK → STRUCTURE
  var plot = document.querySelector('.js-plotly-plot');

  if (plot) {

    plot.on('plotly_click', function(data) {

      if (!data.points || data.points.length === 0) return;

      var site = data.points[0].x;

      console.log('🟡 Clicked site:', site);

      if (!stage) {
        console.log('⚠️ stage not available');
        return;
      }

      var sele = site + ' and (:A or :B)';

      stage.compList.forEach(function(comp) {

        // ✅ highlight clicked residue
        comp.addRepresentation('spacefill', {
          sele: sele,
          color: 'green',
          radius: 2.0
        });

        // ✅ zoom into residue
        comp.autoView(sele);
      });

    });

  }

	// ✅ CONNECT DROPDOWN → STRUCTURE
  var plot = document.querySelector('.js-plotly-plot');

  if (plot) {

plot.on('plotly_buttonclicked', function(ev) {

  if (!ev || !ev.button || !ev.button.label) return;

  var label = ev.button.label;

  console.log('🔄 Dropdown switched to:', label);

  // ✅ ignore undefined / duplicate / invalid
  if (!label || !(label in window.highlight_sites)) {
    console.log('⛔ Ignoring invalid dropdown event');
    return;
  }

  if (window.current_receptor === label) {
    console.log('⛔ Duplicate event ignored');
    return;
  }

  window.current_receptor = label;

  loadStructure(label);
});
  }

});
