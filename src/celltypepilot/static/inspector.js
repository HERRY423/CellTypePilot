// CellTypePilot — Web Inspector frontend logic
// Depends on: `evidenceData` (injected by the dashboard template)

const overrides = {};

function filterTable() {
    const confFilter = document.getElementById('confidenceFilter').value;
    const flagFilter = document.getElementById('flagFilter').value;
    const search = document.getElementById('searchBox').value.toLowerCase();
    const rows = document.querySelectorAll('#annotationTable tbody tr');

    rows.forEach(row => {
        const conf = row.dataset.confidence;
        const flag = row.dataset.flag;
        const ct = row.dataset.celltype;
        const show = (!confFilter || conf === confFilter)
                  && (!flagFilter || flag === flagFilter)
                  && (!search || ct.includes(search));
        row.style.display = show ? '' : 'none';
    });
}

// Helper: build a <p> with a bold label + text content (no innerHTML)
function _makeTextParagraph(label, value) {
    const p = document.createElement('p');
    const strong = document.createElement('strong');
    strong.textContent = label;
    p.appendChild(strong);
    p.appendChild(document.createTextNode(' ' + value));
    return p;
}

function showEvidence(cluster) {
    document.getElementById('evidenceCluster').textContent = cluster;
    const ev = evidenceData[cluster] || {};
    const container = document.getElementById('evidenceContent');
    container.textContent = '';

    container.appendChild(_makeTextParagraph('Cell Type:', ev.cell_type || 'N/A'));
    container.appendChild(_makeTextParagraph('Score:', (ev.combined_score || 0).toFixed(3)));
    container.appendChild(_makeTextParagraph('Marker Overlap:', ((ev.pct_overlap || 0) * 100).toFixed(0) + '%'));

    const h4 = document.createElement('h4');
    h4.style.marginTop = '12px';
    h4.textContent = 'Top Markers';
    container.appendChild(h4);

    const ul = document.createElement('ul');
    ul.className = 'evidence-list';
    (ev.top_markers || []).forEach(m => {
        const li = document.createElement('li');
        const gene = document.createElement('span');
        gene.className = 'gene';
        gene.textContent = m.gene;
        li.appendChild(gene);
        li.appendChild(document.createTextNode(
            ' (pct=' + (m.pct * 100).toFixed(0) + '%, FC=' + m.fc.toFixed(2) + ')'
        ));
        ul.appendChild(li);
    });
    container.appendChild(ul);

    if (ev.critic_notes) {
        const p = _makeTextParagraph('Critic Notes:', ev.critic_notes);
        p.style.marginTop = '12px';
        container.appendChild(p);
    }

    document.getElementById('evidenceModal').style.display = 'block';
}

function showOverride(cluster, currentType) {
    document.getElementById('overrideCluster').textContent = cluster;
    document.getElementById('currentAnnotation').textContent = currentType;
    document.getElementById('newCellType').value = overrides[cluster]?.new_type || '';
    document.getElementById('overrideReason').value = overrides[cluster]?.reason || '';
    document.getElementById('overrideModal').style.display = 'block';
}

function saveOverride() {
    const cluster = document.getElementById('overrideCluster').textContent;
    const newType = document.getElementById('newCellType').value;
    const reason = document.getElementById('overrideReason').value;
    if (!newType) { alert('Please enter a new cell type'); return; }

    // Send to server (not just browser memory)
    fetch('/api/override', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({cluster: cluster, new_type: newType, reason: reason})
    })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            overrides[cluster] = { new_type: newType, reason: reason };
            closeModal('overrideModal');
            showSaveBar();
        } else {
            alert('Error: ' + (data.error || 'Unknown error'));
        }
    })
    .catch(err => alert('Network error: ' + err));
}

function showSaveBar() {
    const bar = document.getElementById('saveBar');
    const count = Object.keys(overrides).length;
    bar.style.display = 'flex';
    document.getElementById('overrideCount').textContent = count;
}

function hideSaveBar() {
    document.getElementById('saveBar').style.display = 'none';
}

function applyOverrides() {
    if (Object.keys(overrides).length === 0) {
        alert('No overrides to apply.');
        return;
    }
    if (!confirm('Apply ' + Object.keys(overrides).length + ' override(s) to .h5ad?\nThis will modify the data file and create a backup.')) return;

    fetch('/api/overrides/apply', { method: 'POST' })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            alert('Applied ' + data.result.applied + ' override(s).\nBackup: ' + data.result.backup);
            location.reload();  // Reload to show updated data
        } else {
            alert('Error: ' + (data.error || 'Apply failed'));
        }
    })
    .catch(err => alert('Network error: ' + err));
}

function exportOverrides() {
    if (Object.keys(overrides).length === 0) {
        alert('No overrides to export.');
        return;
    }
    // Download from server (authoritative copy)
    fetch('/api/overrides')
    .then(r => r.json())
    .then(data => {
        const blob = new Blob([JSON.stringify(data.overrides, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'annotation_overrides.json';
        a.click();
        URL.revokeObjectURL(url);
    });
}

function closeModal(id) {
    document.getElementById(id).style.display = 'none';
}

document.addEventListener('DOMContentLoaded', () => {
    // Close modals on overlay click
    document.querySelectorAll('.modal-overlay').forEach(el => {
        el.addEventListener('click', e => { if (e.target === el) el.style.display = 'none'; });
    });

    // Load existing overrides from server on page load
    fetch('/api/overrides').then(r => r.json()).then(data => {
        Object.assign(overrides, data.overrides || {});
        if (Object.keys(overrides).length > 0) showSaveBar();
    });
});
