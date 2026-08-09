// CellTypePilot — Web Inspector Cockpit frontend logic
// Depends on: `evidenceData` (injected by the dashboard template)

const overrides = {};
let activeClusterForHistory = null;

function switchTab(tabId) {
    document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.nav-tab').forEach(el => el.classList.remove('active'));

    const targetTab = document.getElementById(tabId);
    if (targetTab) targetTab.classList.add('active');

    const navBtn = document.querySelector(`.nav-tab[data-tab="${tabId}"]`);
    if (navBtn) navBtn.classList.add('active');
}

function filterTable() {
    const confFilter = document.getElementById('confidenceFilter')?.value || '';
    const flagFilter = document.getElementById('flagFilter')?.value || '';
    const statusFilter = document.getElementById('statusFilter')?.value || '';
    const search = (document.getElementById('searchBox')?.value || '').toLowerCase();
    const rows = document.querySelectorAll('#annotationTable tbody tr');

    rows.forEach(row => {
        const conf = row.dataset.confidence || '';
        const flag = row.dataset.flag || '';
        const status = row.dataset.status || '';
        const ct = row.dataset.celltype || '';
        const show = (!confFilter || conf === confFilter)
                  && (!flagFilter || flag === flagFilter)
                  && (!statusFilter || status === statusFilter)
                  && (!search || ct.includes(search));
        row.style.display = show ? '' : 'none';
    });
}

function _makeTextParagraph(label, value) {
    const p = document.createElement('p');
    const strong = document.createElement('strong');
    strong.textContent = label;
    p.appendChild(strong);
    p.appendChild(document.createTextNode(' ' + value));
    return p;
}

function _markerList(title, genes, cssClass) {
    const wrap = document.createElement('div');
    wrap.style.marginTop = '10px';
    const h = document.createElement('h4');
    h.textContent = title;
    wrap.appendChild(h);
    if (!genes || !genes.length) {
        const p = document.createElement('p');
        p.className = 'text-muted';
        p.textContent = '—';
        wrap.appendChild(p);
        return wrap;
    }
    const ul = document.createElement('ul');
    ul.className = 'evidence-list ' + (cssClass || '');
    genes.forEach(g => {
        const li = document.createElement('li');
        const span = document.createElement('span');
        span.className = 'gene';
        span.textContent = g;
        li.appendChild(span);
        ul.appendChild(li);
    });
    wrap.appendChild(ul);
    return wrap;
}

function showEvidence(cluster) {
    document.getElementById('evidenceCluster').textContent = cluster;
    const container = document.getElementById('evidenceContent');
    container.textContent = '';
    const loading = document.createElement('p');
    loading.className = 'text-muted';
    loading.textContent = 'Loading Identity × State × Novelty panel…';
    container.appendChild(loading);
    document.getElementById('evidenceModal').style.display = 'block';

    fetch('/api/clusters/' + encodeURIComponent(cluster) + '/review-panel')
    .then(r => r.json())
    .then(panel => {
        container.textContent = '';
        if (!panel || panel.ok === false) {
            // Fallback to embedded evidenceData
            const ev = evidenceData[cluster] || {};
            container.appendChild(_makeTextParagraph('Cell Type:', ev.cell_type || 'N/A'));
            container.appendChild(_makeTextParagraph('Novelty/OOD:', (ev.novelty_decision || 'not_assessed')));
            return;
        }
        const axes = panel.axes || {};
        const id = axes.identity || {};
        const st = axes.state || {};
        const nov = axes.novelty || {};

        const grid = document.createElement('div');
        grid.className = 'review-axes-grid';
        [
            ['Identity', id, [
                ['Label', id.cell_type],
                ['Decision', id.decision],
                ['CL ID', id.cl_id],
                ['Score', (id.evidence_score || 0).toFixed(3)],
                ['Critic', id.critic_confidence + ' / ' + id.critic_flags],
            ]],
            ['State', st, [
                ['Candidate', st.state_candidate],
                ['Decision', st.state_decision],
                ['Score', (st.state_score || 0).toFixed(3)],
                ['Confidence', st.state_confidence],
            ]],
            ['Novelty / OOD', nov, [
                ['Decision', nov.novelty_decision],
                ['Score', (nov.novelty_score || 0).toFixed(3)],
                ['Unmapped', (nov.top_unmapped_markers || []).join(', ') || '—'],
            ]],
        ].forEach(([title, _axis, rows]) => {
            const card = document.createElement('div');
            card.className = 'review-axis-card';
            const h = document.createElement('h4');
            h.textContent = title;
            card.appendChild(h);
            rows.forEach(([k, v]) => card.appendChild(_makeTextParagraph(k + ':', v == null || v === '' ? '—' : String(v))));
            grid.appendChild(card);
        });
        container.appendChild(grid);

        container.appendChild(_markerList('Supporting markers', id.supporting_markers, 'markers-support'));
        container.appendChild(_markerList('Opposing / negative-expressed', id.opposing_markers, 'markers-oppose'));
        container.appendChild(_markerList('Silent expected markers', id.silent_markers, 'markers-silent'));
        container.appendChild(_markerList('Missing expected markers', id.missing_markers, 'markers-missing'));

        const neigh = id.neighbor_candidates || [];
        if (neigh.length) {
            const h = document.createElement('h4');
            h.style.marginTop = '12px';
            h.textContent = 'Neighbor / runner-up candidates';
            container.appendChild(h);
            const ul = document.createElement('ul');
            ul.className = 'evidence-list';
            neigh.forEach(n => {
                const li = document.createElement('li');
                li.textContent = (n.cell_type || '') + (n.score != null ? ' (score=' + Number(n.score).toFixed(3) + ')' : '');
                ul.appendChild(li);
            });
            container.appendChild(ul);
        }

        const strata = panel.donor_batch_strata || {};
        const sh = document.createElement('h4');
        sh.style.marginTop = '12px';
        sh.textContent = 'Donor / batch / sample strata';
        container.appendChild(sh);
        ['donors', 'batches', 'samples'].forEach(key => {
            const block = strata[key] || {};
            const status = block.status || 'not_assessed';
            const levels = (block.levels || []).map(l => l.label + '=' + l.n_cells).join(', ') || '—';
            container.appendChild(_makeTextParagraph(key + ' (' + status + '):', levels));
        });

        const lit = panel.literature || {};
        container.appendChild(_makeTextParagraph(
            'Literature / provenance:',
            (lit.status || 'not_assessed') +
            (lit.pmids && lit.pmids.length ? ' PMIDs: ' + lit.pmids.join(', ') : '') +
            (lit.sources && lit.sources.length ? ' sources: ' + lit.sources.join(', ') : '')
        ));
        if (lit.note) {
            const n = document.createElement('p');
            n.className = 'text-muted';
            n.textContent = lit.note;
            container.appendChild(n);
        }

        if (id.critic_notes) {
            container.appendChild(_makeTextParagraph('Critic notes:', id.critic_notes));
        }
        if (nov.alternative_explanations) {
            container.appendChild(_makeTextParagraph('Alternatives:', nov.alternative_explanations));
        }
        if (nov.recommended_next_actions) {
            container.appendChild(_makeTextParagraph('Next actions:', nov.recommended_next_actions));
        }

        const policy = panel.edit_policy || {};
        const foot = document.createElement('p');
        foot.className = 'text-muted';
        foot.style.marginTop = '12px';
        foot.textContent = 'Edits are append-only; after Apply, derived artifacts go stale until regenerate + resign. ' +
            'append_only=' + !!policy.append_only_audit +
            ', stale_after_apply=' + !!policy.derived_artifacts_stale_after_apply;
        container.appendChild(foot);
    })
    .catch(err => {
        container.textContent = '';
        container.appendChild(_makeTextParagraph('Error loading panel:', String(err)));
    });
}

function resignReview() {
    const signer = prompt('Signer identity for re-sign', 'web_reviewer');
    if (!signer) return;
    fetch('/api/review/resign', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({signer: signer, regenerate: true})
    })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            alert('Re-signed. Artifacts current. signature=' + (data.signature?.signature_sha256 || '').slice(0, 16) + '…');
            location.reload();
        } else {
            alert('Resign failed: ' + (data.error || 'unknown'));
        }
    })
    .catch(err => alert('Network error: ' + err));
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
            location.reload();
        } else {
            alert('Error: ' + (data.error || 'Unknown error'));
        }
    })
    .catch(err => alert('Network error: ' + err));
}

function updateClusterStatus(clusterId, status) {
    fetch(`/api/clusters/${clusterId}/status`, {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({status: status})
    })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            const row = document.querySelector(`#annotationTable tbody tr[data-celltype][data-confidence]`);
            // Update checklist & audit trail
            refreshAudit();
            fetchChecklist();
        } else {
            alert('Status update failed: ' + (data.error || 'Unknown error'));
        }
    });
}

function showClusterHistory(clusterId) {
    activeClusterForHistory = clusterId;
    document.getElementById('historyClusterId').textContent = clusterId;
    document.getElementById('historyModal').style.display = 'block';

    fetch(`/api/clusters/${clusterId}/history`)
    .then(r => r.json())
    .then(data => {
        if (!data.ok) return;

        // Baseline Card
        const base = data.baseline || {};
        const baseBox = document.getElementById('historyBaseline');
        baseBox.textContent = '';
        baseBox.appendChild(_makeTextParagraph('Original Model Assignment:', base.cell_type || 'Unknown'));
        baseBox.appendChild(_makeTextParagraph('CL ID:', base.cl_id || '-'));
        baseBox.appendChild(_makeTextParagraph('Evidence Score:', (base.combined_score || 0).toFixed(3)));
        baseBox.appendChild(_makeTextParagraph('Critic Flags:', base.critic_flags || 'PASS'));

        // Notes List
        const notesBox = document.getElementById('historyNotesList');
        notesBox.textContent = '';
        const notes = data.notes || [];
        if (notes.length === 0) {
            notesBox.textContent = 'No review notes recorded yet.';
        } else {
            notes.forEach(n => {
                const item = document.createElement('div');
                item.className = 'note-item';

                const author = document.createElement('span');
                author.className = 'note-author';
                author.textContent = n.author || 'Reviewer';
                item.appendChild(author);

                const time = document.createElement('span');
                time.className = 'note-time';
                time.textContent = (n.timestamp || '').substring(0, 19);
                item.appendChild(time);

                const textP = document.createElement('p');
                textP.textContent = n.text;
                textP.style.marginTop = '4px';
                item.appendChild(textP);

                notesBox.appendChild(item);
            });
        }

        // Audit Trail for cluster
        const auditList = document.getElementById('historyAuditList');
        auditList.textContent = '';
        const events = data.audit_history || [];
        if (events.length === 0) {
            const li = document.createElement('li');
            li.textContent = 'No cluster events recorded';
            auditList.appendChild(li);
        } else {
            events.forEach(ev => {
                const li = document.createElement('li');
                li.textContent = `${(ev.timestamp || '').substring(0, 19)} — ${ev.event_type}`;
                auditList.appendChild(li);
            });
        }
    });
}

function addClusterNote() {
    if (!activeClusterForHistory) return;
    const author = document.getElementById('noteAuthor').value.trim() || 'Reviewer';
    const text = document.getElementById('noteText').value.trim();

    if (!text) { alert('Please enter note text'); return; }

    fetch(`/api/clusters/${activeClusterForHistory}/note`, {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({author: author, text: text})
    })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            document.getElementById('noteText').value = '';
            showClusterHistory(activeClusterForHistory);
        } else {
            alert('Add note failed: ' + (data.error || 'Unknown error'));
        }
    });
}

function toggleChecklistItem(itemKey, completed) {
    fetch('/api/checklist', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({item_key: itemKey, completed: completed})
    })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            updateChecklistUI(data.checklist);
        }
    });
}

function fetchChecklist() {
    fetch('/api/checklist')
    .then(r => r.json())
    .then(data => {
        if (data.ok) updateChecklistUI(data.checklist);
    });
}

function updateChecklistUI(chk) {
    if (!chk) return;
    const pct = chk.readiness_pct || 0;
    const bar = document.getElementById('checklistProgressBar');
    if (bar) bar.style.width = pct + '%';

    const pctText = document.getElementById('checklistPct');
    if (pctText) pctText.textContent = pct + '%';

    const pill = document.getElementById('checklistTabPill');
    if (pill) pill.textContent = pct + '%';
}

function submitSignoff(e) {
    e.preventDefault();
    const name = document.getElementById('reviewerName').value.trim();
    const role = document.getElementById('reviewerRole').value.trim();
    const decision = document.getElementById('reviewDecision').value;
    const notes = document.getElementById('signoffNotes').value.trim();
    const force = document.getElementById('forceSignoff').checked;

    fetch('/api/signoff', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            reviewer_name: name,
            reviewer_role: role,
            decision: decision,
            notes: notes,
            force: force
        })
    })
    .then(r => r.json())
    .then(data => {
        if (data.ok) {
            alert('Sign-off certificate successfully created!');
            location.reload();
        } else {
            alert('Sign-off failed: ' + (data.error || 'Unknown error'));
        }
    });
}

function showSaveBar() {
    const bar = document.getElementById('saveBar');
    const count = Object.keys(overrides).length;
    if (bar) {
        bar.style.display = 'flex';
        document.getElementById('overrideCount').textContent = count;
    }
}

function hideSaveBar() {
    const bar = document.getElementById('saveBar');
    if (bar) bar.style.display = 'none';
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
            location.reload();
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

function _renderList(id, values, emptyText) {
    const list = document.getElementById(id);
    if (!list) return;
    list.textContent = '';
    if (!values || values.length === 0) {
        const li = document.createElement('li');
        li.textContent = emptyText;
        list.appendChild(li);
        return;
    }
    values.forEach(value => {
        const li = document.createElement('li');
        li.textContent = value;
        list.appendChild(li);
    });
}

function refreshAudit() {
    fetch('/api/artifact-status')
    .then(r => r.json())
    .then(data => {
        const status = data.artifact_status || {};
        const state = document.getElementById('artifactState');
        const message = document.getElementById('artifactMessage');
        if (state) state.textContent = status.review_state || 'unknown';
        if (message) message.textContent = status.message || '';
        _renderList('staleArtifacts', status.stale_artifacts || [], 'None');
    });

    fetch('/api/audit?limit=5')
    .then(r => r.json())
    .then(data => {
        const events = (data.events || []).map(ev => (ev.timestamp || '').substring(0, 19) + ' — ' + (ev.event_type || 'unknown'));
        _renderList('auditEvents', events, 'No review events yet');
    });
}

function closeModal(id) {
    document.getElementById(id).style.display = 'none';
}

// ── Read-only run observability (never mutates predictions) ──────────
let _obsTimer = null;

function toggleObsAutoRefresh(enabled) {
    if (_obsTimer) {
        clearInterval(_obsTimer);
        _obsTimer = null;
    }
    if (enabled) {
        _obsTimer = setInterval(refreshObservability, 10000);
    }
}

function _fmtNum(value, digits) {
    if (value === null || value === undefined || Number.isNaN(value)) return '—';
    return Number(value).toFixed(digits);
}

function refreshObservability() {
    fetch('/api/observability')
    .then(r => r.json())
    .then(data => {
        const obs = data.observability || {};
        if (obs.prediction_mutation_allowed) {
            console.warn('Server unexpectedly allowed prediction mutation on observability surface');
        }
        const cp = obs.checkpoints || {};
        const eta = obs.fold_eta || {};
        const host = obs.host || {};
        const cpu = host.cpu || {};
        const gpu = host.gpu || {};
        const stale = obs.stale || {};

        const nEl = document.getElementById('obsNCheckpoints');
        if (nEl) nEl.textContent = cp.n_status_files || 0;
        const prog = document.getElementById('obsProgress');
        if (prog) {
            prog.textContent = `${eta.n_completed || 0} / ${eta.n_running || 0} / ${eta.n_failed || 0}`;
        }
        const etaEl = document.getElementById('obsEta');
        if (etaEl) {
            etaEl.textContent = eta.estimated_remaining_seconds == null
                ? '—'
                : _fmtNum(eta.estimated_remaining_seconds, 0);
        }
        const staleEl = document.getElementById('obsStale');
        if (staleEl) staleEl.textContent = stale.derived_artifacts_stale ? 'yes' : 'no';

        const rootEl = document.getElementById('obsRunRoot');
        if (rootEl) rootEl.textContent = obs.run_root || '';

        const hostBody = document.getElementById('obsHostBody');
        if (hostBody) {
            hostBody.textContent = '';
            hostBody.appendChild(_makeTextParagraph('CPU%:', String(cpu.cpu_percent ?? 'n/a') + ' · cores ' + String(cpu.cpu_count_logical ?? 'n/a')));
            hostBody.appendChild(_makeTextParagraph('Platform:', String(cpu.platform || 'n/a')));
            if (gpu.available && (gpu.devices || []).length) {
                hostBody.appendChild(_makeTextParagraph('GPU source:', String(gpu.source || '')));
                const ul = document.createElement('ul');
                ul.className = 'audit-list';
                (gpu.devices || []).forEach(dev => {
                    const li = document.createElement('li');
                    let text = `${dev.index}: ${dev.name || ''}`;
                    if (dev.utilization_gpu_percent != null) text += ` · util ${dev.utilization_gpu_percent}%`;
                    if (dev.memory_used_mib != null) text += ` · mem ${dev.memory_used_mib}/${dev.memory_total_mib} MiB`;
                    li.textContent = text;
                    ul.appendChild(li);
                });
                hostBody.appendChild(ul);
            } else {
                const p = document.createElement('p');
                p.className = 'text-muted';
                p.textContent = 'No GPU reported (nvidia-smi / torch.cuda unavailable).';
                hostBody.appendChild(p);
            }
        }

        const staleBody = document.getElementById('obsStaleBody');
        if (staleBody) {
            staleBody.textContent = '';
            staleBody.appendChild(_makeTextParagraph('Source:', String(stale.source || 'n/a')));
            staleBody.appendChild(_makeTextParagraph('Review state:', String(stale.review_state || 'n/a')));
            staleBody.appendChild(_makeTextParagraph('Reason:', String(stale.stale_reason || 'n/a')));
            const p = document.createElement('p');
            p.className = 'text-muted';
            p.appendChild(document.createTextNode('Run root: '));
            const code = document.createElement('code');
            code.id = 'obsRunRoot';
            code.textContent = obs.run_root || '';
            p.appendChild(code);
            staleBody.appendChild(p);
        }

        const tbody = document.getElementById('obsCheckpointBody');
        if (tbody) {
            tbody.textContent = '';
            const rows = (cp.records || []);
            if (!rows.length) {
                const tr = document.createElement('tr');
                const td = document.createElement('td');
                td.colSpan = 6;
                td.className = 'text-muted';
                td.textContent = 'No checkpoint status files found under this run root.';
                tr.appendChild(td);
                tbody.appendChild(tr);
            } else {
                rows.forEach(row => {
                    const tr = document.createElement('tr');
                    tr.className = 'obs-status-' + (row.status || 'unknown');
                    const cells = [
                        row.method,
                        row.fold_id,
                        row.status,
                        `${_fmtNum(row.elapsed_seconds, 1)} / ${_fmtNum(row.duration_seconds, 1)}`,
                        row.failure_reason || row.error || '—',
                        row.prediction_sha256 ? (row.prediction_sha256.slice(0, 16) + '…') : '—',
                    ];
                    cells.forEach((value, idx) => {
                        const td = document.createElement('td');
                        if (idx === 1 || idx === 5) {
                            const code = document.createElement('code');
                            if (idx === 5) code.className = 'hash';
                            code.textContent = value;
                            td.appendChild(code);
                        } else {
                            td.textContent = value;
                        }
                        tr.appendChild(td);
                    });
                    tbody.appendChild(tr);
                });
            }
        }

        const failList = document.getElementById('obsFailures');
        if (failList) {
            const failures = obs.failures || [];
            _renderList(
                'obsFailures',
                failures.map(f => `${f.method} / ${f.fold_id}: ${f.failure_reason}`),
                'No failed checkpoints.',
            );
        }

        const pbody = document.getElementById('obsProductBody');
        if (pbody) {
            pbody.textContent = '';
            const products = obs.products || [];
            if (!products.length) {
                const tr = document.createElement('tr');
                const td = document.createElement('td');
                td.colSpan = 4;
                td.className = 'text-muted';
                td.textContent = 'No product inventory.';
                tr.appendChild(td);
                pbody.appendChild(tr);
            } else {
                products.forEach(p => {
                    const tr = document.createElement('tr');
                    [p.name, p.present ? 'yes' : 'no', p.byte_size ?? '—', p.sha256 || p.hash_skipped || p.error || '—']
                        .forEach((value, idx) => {
                            const td = document.createElement('td');
                            if (idx === 0 || idx === 3) {
                                const code = document.createElement('code');
                                if (idx === 3) code.className = 'hash';
                                code.textContent = value;
                                td.appendChild(code);
                            } else {
                                td.textContent = value;
                            }
                            tr.appendChild(td);
                        });
                    pbody.appendChild(tr);
                });
            }
        }

        const pill = document.getElementById('obsTabPill');
        if (pill) pill.textContent = String(cp.n_status_files || 0);
    })
    .catch(err => console.warn('observability refresh failed', err));
}

document.addEventListener('DOMContentLoaded', () => {
    document.querySelectorAll('.modal-overlay').forEach(el => {
        el.addEventListener('click', e => { if (e.target === el) el.style.display = 'none'; });
    });

    fetch('/api/overrides').then(r => r.json()).then(data => {
        Object.assign(overrides, data.overrides || {});
        if (Object.keys(overrides).length > 0) showSaveBar();
    });

    refreshAudit();
    fetchChecklist();
    // Initial observability hydrate (read-only)
    refreshObservability();
});
