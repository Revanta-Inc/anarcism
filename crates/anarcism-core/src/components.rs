pub(crate) fn connected_components<T>(
    items: Vec<T>,
    linked: impl Fn(&T, &T) -> bool,
) -> Vec<Vec<T>> {
    let mut parents: Vec<usize> = (0..items.len()).collect();
    for left in 0..items.len() {
        for right in left + 1..items.len() {
            if linked(&items[left], &items[right]) {
                union(&mut parents, left, right);
            }
        }
    }

    let roots: Vec<_> = (0..items.len())
        .map(|index| find(&mut parents, index))
        .collect();
    let mut group_by_root = vec![usize::MAX; items.len()];
    let mut groups = Vec::new();
    for (item, root) in items.into_iter().zip(roots) {
        let group_index = if group_by_root[root] == usize::MAX {
            group_by_root[root] = groups.len();
            groups.push(Vec::new());
            groups.len() - 1
        } else {
            group_by_root[root]
        };
        groups[group_index].push(item);
    }
    groups
}

fn find(parents: &mut [usize], index: usize) -> usize {
    if parents[index] != index {
        parents[index] = find(parents, parents[index]);
    }
    parents[index]
}

fn union(parents: &mut [usize], left: usize, right: usize) {
    let left_root = find(parents, left);
    let right_root = find(parents, right);
    if left_root != right_root {
        parents[right_root] = left_root;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bridging_item_merges_components() {
        let groups = connected_components(
            vec![(0_usize, 20_usize), (14, 100), (15, 25)],
            |left, right| {
                let overlap = left.1.min(right.1).saturating_sub(left.0.max(right.0));
                let shorter = (left.1 - left.0).min(right.1 - right.0);
                overlap * 2 >= shorter
            },
        );
        assert_eq!(groups.len(), 1);
    }
}
